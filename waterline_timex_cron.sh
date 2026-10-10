#!/bin/bash
# Timex waterline capture, archive, and mapping -- safe to run repeatedly
# ------------------------------------------------------------------
# Runs the detector over whatever timex images are currently present,
# archives the results permanently, matches them to measured water
# level, and regenerates the elevation maps.
#
# WHY IT IS BUILT THIS WAY:
#   * waterline_detector_v5.py CLEARS its output folder at the start of
#     every run, so the working folder can never be the accumulator --
#     each run would destroy the previous day's detections. The
#     detector writes to a scratch folder and results are copied into
#     archive/processed_timex with `cp -n` (no clobber), which
#     accumulates and is idempotent.
#   * Because it is idempotent, running it more often than necessary is
#     harmless. That matters: cleanup.sh runs at 19:45 AND 21:00 local
#     while captures continue to 20:59, so a single nightly run cannot
#     cover the day.
#   * cleanup.sh uploads to S3 BEFORE deleting, and only deletes after
#     the upload is confirmed, so a missed run means "recover from S3",
#     not data loss.
#
# Crontab entries (LOCAL time, matching the rest of that file):
#   30 12 * * * /mnt/I2Rgus_Data/waterline/waterline_timex_cron.sh
#   25 19 * * * /mnt/I2Rgus_Data/waterline/waterline_timex_cron.sh
#   55 20 * * * /mnt/I2Rgus_Data/waterline/waterline_timex_cron.sh

set -u

BASE=/mnt/I2Rgus_Data/waterline
SCRATCH="$BASE/processed_timex"
# The detector also CLEARS its input folder every run. Its default,
# $BASE/input, is where collect_all_ground_truth.py looks for images to
# click, so this job must never use it -- it would wipe images staged
# for ground-truth collection three times a day.
SCRATCH_INPUT="$BASE/input_timex"
LOCK="$BASE/.timex_cron.lock"
ARCHIVE_CSV="$BASE/archive/processed_timex"
ARCHIVE_IMG="$BASE/archive/images_timex"
ARCHIVE_RAS="$BASE/archive/ras_c2"
LOG="$BASE/logs/timex_cron.log"

# Keep every timex image. Disk is not the constraint (1.7 TB free, this
# grows ~9 GB/year) and keeping them is what allows the whole archive
# to be REPROCESSED when the detector or its thresholds change -- which
# has already been necessary several times. Set to a glob such as
# "*12_00_00*.timex.jpg" to keep only a daily backdrop instead.
BACKDROP_PATTERN="all"

# The other Argus products, kept the same way, each in archive/images_<product>.
# The wave models read the BRIGHT image, so without it the live wave record
# can never be recomputed when a model improves; dark/snap/var are what a
# multi-product model would train on. About 4.5 GB/year per product per camera.
# Set to "" to stop (what is already archived stays).
ARCHIVE_PRODUCTS="bright dark snap var"

# Water level source. GNSS-R is preferred over the tide model: it is a
# MEASUREMENT at this station, already referenced to NAVD88 (so it
# needs no datum offset), and it includes surge and setup that an
# astronomical model cannot predict. Its cost is latency -- gnssrefl
# processing runs about 2 days behind, so the newest maps cover
# day-2 through day-8 rather than today. daily_elevation_map.py
# anchors its window to the last date present in the data, not to
# today, so that lag is reported honestly rather than silently
# shortening the range.
GNSSR_SPLINE=/home/argus_user/GNSS/v4.1/products/refl_code/Files/usgs/usgs_spline_out.txt
CONTOURS="$BASE/contour_points_timex.csv"

# Waterline consistency (waterline_consistency.py). In one image column
# the waterline row must move DOWN the photo as the water rises (nearer
# the camera); lines, or parts of lines, that the other lines from the
# same camera and nearby days clearly contradict (by more than 0.75 m of
# elevation) are left out of the maps, the ground points and the DEM. On
# 29 Sep - 5 Oct 2026 C1 drew +1 m lines seaward of -0.5 m lines on the
# right of its frame, and the DEM's southern ~80 m (C1 only) came out
# blanked for spread. Columns where the rows run AGAINST the water level
# -- the detector wrong for most lines there, e.g. a search-envelope
# floor above the high-tide waterline -- are named in a WARNING in this
# log, and the lines there are judged against the lower-water lines.
# The unfiltered file is kept as it is; the filtered copy is
# CONTOURS_QC, the report of every line affected is CONSISTENCY_REPORT,
# and waterline_consistency_<cam>.png shows the last week's rejections on
# a photo. If the filter fails, everything downstream uses the unfiltered
# file (WARNING in the log) -- no product is lost. Set CONSISTENCY_ENABLE=0
# to skip it.
CONSISTENCY_ENABLE=1
CONTOURS_QC="$BASE/contour_points_timex_qc.csv"
CONSISTENCY_REPORT="$BASE/waterline_consistency_report.csv"
CONSISTENCY_PLOT="$BASE/waterline_consistency"

# Map windows. The 7-day map is the diagnostic one: over a week a
# filled band is still mostly measurement repeatability. Over 30 days
# real morphological change dominates, so the monthly map shows change
# but no longer says anything useful about repeatability.
#
# MAP_MONTH_MAX_LINES caps how many lines the monthly map draws. A
# 30-day window holds 500+ waterlines, which saturate the intertidal
# zone into a solid mass with the photo invisible behind it. The cap
# subsamples EVENLY IN TIME, not every Nth file, because captures are
# irregular -- nights and bad-weather frames are missing, so
# index-based thinning would over-represent the good stretches.
MAP_WEEK_DAYS=7
MAP_MONTH_DAYS=30
MAP_MONTH_MAX_LINES=0

# Georectification and DEM.
#
# The DEM is rebuilt from the FULL archive every run, not just the
# newest points: repeat tidal crossings are what make each cell well
# determined, so more history means a better surface, not a stale one.
#
# DEM_CELL of 2 m is a compromise. Below the georectification's own
# 0.33 m accuracy there is nothing to gain; well above it the beach
# profile gets smoothed away. DEM_MIN_POINTS rejects cells resting on
# one or two detections, and DEM_MAX_SPREAD blanks cells where repeat
# crossings disagree by more than half a metre -- those are not
# measuring a single surface.
#
# Set DEM_ENABLE=0 to skip this stage without editing anything else.
DEM_ENABLE=1
CAL="$BASE/calibration"
GROUND="$BASE/contour_points_ground.csv"
DEM_STEM="$BASE/dem_intertidal"
DEM_CELL=2.0
DEM_MIN_POINTS=3
DEM_MAX_SPREAD=0.5

# Runup from C2 timestacks (runup_from_timestack.py): swash positions
# along each line, in metres.
#
# cleanup.sh moves ras.tiff files to S3, but GNSS-R lags ~2 days, so
# the C2 stacks are copied into $ARCHIVE_RAS (~9 MB each, ~270 MB/day,
# ~100 GB/year) and the runup script runs over that archive with
# --require-water-level: each stack is processed once, on the first run
# after GNSS-R covers it. C2 lines 1 and 2 cross the swash; 3 and 4 sit
# on the upper beach and dune foot. Results: runup_c2_line<N>.csv.
#
# Set RUNUP_ENABLE=0 to skip this stage (stacks are still archived).
RUNUP_ENABLE=1
RUNUP_LINES="1 2"
RAS_SOURCE=/mnt/I2Rgus_Data/ImageProducts/products

# Offshore waves (fetch_buoy_waves.py). NDBC keeps only 45 days online,
# so every run merges the buoy's recent record into WAVES_CSV. Each
# waterline is then tagged with the offshore wave height at its capture
# time, and the DEM leaves out frames above DEM_MAX_HS: in big waves the
# swash sits above still water by the wave setup, so those waterlines
# are given an elevation that is too low for where they lie. On Sep 23
# 2026 frames in 3.8-3.9 m waves (week median 1.0 m) plotted ~0.5 m high.
# 44008 faces the open Atlantic like Marconi; heights are an offshore
# index, not the breaking height. Set DEM_MAX_HS="" to keep all frames.
# GNSS-R quality control (gnssr_qc.py). GNSS-R degrades in big waves
# (Sep 25-26 2026: up to ~4 m above the tide). Each run archives the NOAA
# Chatham gauge in NAVD88, checks every GNSS-R reading against it
# (plus gross-range, spike and rate tests), writes a report and a 7-day
# plot, and the contour and runup stages drop failed readings.
TIDE_GAUGE=8447435
GAUGE_CSV="$BASE/archive/gauge_${TIDE_GAUGE}.csv"
GNSSR_QC_CSV="$BASE/archive/gnssr_qc.csv"

WAVE_BUOY=44008
WAVES_CSV="$BASE/archive/waves_${WAVE_BUOY}.csv"
DEM_MAX_HS=1.5

# Day consistency: leave out a camera-day whose waterlines sit, on median,
# more than this far (m) from the DEM built without that day. Catches days
# whose lines are displaced as a whole, which the tide-direction check
# (direction only) cannot see. Set to "" to disable.
DEM_MAX_DAY_OFFSET=0.15

# Rolling-window DEMs. The whole-archive DEM pools weeks of data, which
# blurs real beach change into "spread" (Sep 2026: 0.278 m for Sep 11-23
# against 0.231 m for Sep 11-14 alone, with more cells filled). Each run
# also builds a DEM from the last DEM_WINDOW_DAYS days, keeps a dated copy
# in DEM_SERIES, and maps the change from the DEM a window earlier
# (dem_change.py). Set DEM_WINDOW_DAYS="" to skip.
DEM_WINDOW_DAYS=7
DEM_SERIES="$BASE/archive/dems"

# Wave-setup correction (extract_elevation_contours.py --setup-coef): each
# waterline's elevation becomes water level + C*sqrt(Hs*L0), Hs and Tp from
# the wave record the contours carry (offshore_hs_m, offshore_tp_s). A timex
# shows the water's edge where waves run up to, above still water, so
# without it the DEM reads low in rough water.
#
# C = 0.037, kept after a recheck in Oct 2026. What each fit on the
# 29 Sep - 5 Oct 2026 waterlines gives:
#  - the first fit (8 Oct), against the 2026-09-29 RTK transects: per frame
#    C = -(line - RTK)/sqrt(Hs*L0), each line taken where it lay WITHOUT
#    setup: median 0.037, 90% 0.028-0.058 (14 frames; every crossing,
#    so the survey-floor bias below applies to it too);
#  - repeat crossings of a cell, no survey: dem_from_contours.py --fit-setup
#    on lines built with C = 0 prints a near-best band (the C whose median
#    spread is within 5 mm of its best) of +0.006 to +0.034 on the filtered
#    lines with --max-hs 1.5 (best 0.026, a flat minimum; 0.037 spreads
#    6 mm more) and +0.026 to +0.040 on the unfiltered rows, also with
#    --max-hs 1.5 (best 0.030); the same spread with each line re-projected
#    at still water + setup gives 0.029 (90% ~0.00-0.04); the within-cell
#    slope of elevation on sqrt(Hs*L0), each line re-projected, gives 0.043
#    (90% 0.026-0.053);
#  - the RTK transects stop at their lowest shot, +1.17 to +1.80 m. A line
#    below that can only read low on them, and one that would read high at
#    that level lands seaward of the survey and drops out, so a fit over
#    every crossing is biased by where the survey stops (re-projected, all
#    crossings: 0.067). Over only lines 0.1-0.2 m or more above the lowest
#    shots, re-projected: ~0.043-0.046 (day-block 90% ~0.02-0.055).
# (The re-projected fits, 0.029, 0.043, 0.067 and 0.043-0.046, come from a
# scratch analysis made in Oct 2026, not from code in this repository.)
# 0.037 lies inside each 90% interval above that the survey floor does not
# bias (~0.00-0.04, 0.026-0.053, ~0.02-0.055) and inside the --fit-setup
# band of the unfiltered rows; the filtered lines' band stops just below
# it. It is the Stockdon et al. (2006) setup, 0.35*beta_f*sqrt(H0*L0), for
# a ~1:10 beach face (beta_f ~0.106). The RTK cannot pin C more closely:
# an RTK survey reaching the low-tide line, on days with different waves,
# would. Refit (repeat crossings on contours built with SETUP_COEF="") if
# the wave source (WAVE_BUOY / USE_MARCONI_WAVES) changes, and read what
# --fit-setup prints: it offers a C to apply only when C = 0 and the C in
# use (this SETUP_COEF) both lie outside its near-best band and the
# minimum is not flat. One week's curve is noisy and can sit near 0 (the
# C = 0 lines of 29 Sep - 2 Oct 2026 alone, filtered, --max-hs 1.5: band
# -0.020 to +0.010, best -0.004, "not distinguishable from C = 0"), so
# adopt a C only when it holds over several weeks of different waves.
# About +0.37 m at Hs 1 m, Tp 8 s; the DEM moves ~0.06 m per 0.01 of C.
# Changing C changes every DEM, but the
# week-to-week change rebuilds both weeks with the current setting
# (dem_change.py --rebuild), so a new C does not show there as change.
# survey_products.py and historical_forcing.py carry the same C (their
# SETUP_COEF): change all three together.
SETUP_COEF="0.037"
# If GNSS-R falls further behind than this, something has stopped --
# 2 days is normal, so this allows generous margin before complaining.
GNSSR_STALE_DAYS=5
# The date the GNSS-R record starts and the MOST days with readings it has
# had, remembered from the runs that accepted it ("YYYY-MM-DD N"). Contours,
# ground points, DEM and runup are all rebuilt from the WHOLE archive and
# drop every frame the spline does not cover, so a spline that suddenly
# starts later (gnssrefl refitted only part of the record -- what a
# calendar-year fit did on 2 January) or has lost days in the middle would
# silently remove those frames. Such a spline is refused and the previous
# products kept. N is a high-water mark, so a few days lost run after run
# add up and are refused too. Delete this file if the record was shortened
# on purpose (it is written again, from that spline, on the next run).
GNSSR_START_FILE="$BASE/archive/gnssr_record_start.txt"
GNSSR_MAX_LOST_DAYS=3

mkdir -p "$SCRATCH" "$SCRATCH_INPUT" "$ARCHIVE_CSV" "$ARCHIVE_IMG" "$ARCHIVE_RAS" "$(dirname "$LOG")"

log() { echo "[$(date -u '+%Y-%m-%dT%H:%M:%SZ')] $*" >> "$LOG"; }

# One run at a time. The detector clears its scratch folders on start,
# so an overlapping run (a slow run meeting the next cron slot, or a
# manual run) would delete the files the first run is still working on.
exec 9>"$LOCK"
if ! flock -n 9; then
    log "SKIPPED: another run holds $LOCK"
    exit 0
fi

log "=== run start ==="

cd "$BASE" || { log "FATAL: cannot cd to $BASE"; exit 1; }

# 1. Detect. Bias correction is ON: timex-specific corrections were
#    derived in v6.8 and the detector selects them by image suffix, so
#    a timex run no longer picks up the snap-tuned values.
python3 waterline_detector_v5.py \
    --image-suffix timex.jpg \
    --input-dir "$SCRATCH_INPUT" \
    --output-dir "$SCRATCH" \
    --debug-dir "$BASE/debug" >> "$LOG" 2>&1
status=$?
if [ $status -ne 0 ]; then
    log "WARNING: detector exited $status (continuing to archive whatever it produced)"
fi

# 2. Accumulate. -n never overwrites, so re-running is safe and a
#    filename already archived is left alone.
before=$(find "$ARCHIVE_CSV" -name '*.csv' 2>/dev/null | wc -l)
cp -pn "$SCRATCH"/*.csv "$ARCHIVE_CSV"/ 2>/dev/null
after=$(find "$ARCHIVE_CSV" -name '*.csv' 2>/dev/null | wc -l)
log "detections archived: +$((after - before)) (total $after)"

# 3. Keep the source images before cleanup.sh moves them to S3.
if [ "$BACKDROP_PATTERN" = "all" ]; then
    cp -pn /mnt/I2Rgus_Data/ImageProducts/*.timex.jpg "$ARCHIVE_IMG"/ 2>/dev/null
else
    cp -pn /mnt/I2Rgus_Data/ImageProducts/$BACKDROP_PATTERN "$ARCHIVE_IMG"/ 2>/dev/null
fi
img_count=$(find "$ARCHIVE_IMG" -name '*.jpg' 2>/dev/null | wc -l)
log "images in archive: $img_count"
for prod in $ARCHIVE_PRODUCTS; do
    mkdir -p "$BASE/archive/images_$prod"
    for src in /mnt/I2Rgus_Data/ImageProducts /mnt/I2Rgus_Data/ImageProducts/products; do
        cp -pn "$src"/*."$prod".jpg "$BASE/archive/images_$prod"/ 2>/dev/null
    done
    log "$prod images in archive: $(find "$BASE/archive/images_$prod" -name '*.jpg' 2>/dev/null | wc -l)"
done

cp -pn "$RAS_SOURCE"/*.c2.ras.tiff "$ARCHIVE_RAS"/ 2>/dev/null
ras_count=$(find "$ARCHIVE_RAS" -name '*.ras.tiff' 2>/dev/null | wc -l)
log "C2 timestacks in archive: $ras_count"

# 4. Flag if today produced nothing -- this is the failure mode that
#    would otherwise go unnoticed, since the maps still render from
#    older archived data and look fine.
# LOCAL date: find -newermt reads the string as local midnight, so a UTC
# date put the reference in the future at the 20:55 run (already the
# next day in UTC) and the warning fired every evening.
today=$(date '+%Y-%m-%d')
recent=$(find "$ARCHIVE_CSV" -name '*.csv' -newermt "$today" 2>/dev/null | wc -l)
if [ "$recent" -eq 0 ]; then
    log "WARNING: no detections archived with today's date ($today). If this repeats,"
    log "         the capture chain or the detector is failing and the maps will"
    log "         quietly go stale rather than erroring."
fi

# 5a. Offshore waves. A failed download is only a warning: the archive
#     still covers everything up to the last successful run.
python3 "$BASE/fetch_buoy_waves.py" --station "$WAVE_BUOY" --output "$WAVES_CSV" >> "$LOG" 2>&1 \
    || log "WARNING: buoy $WAVE_BUOY download failed; using the existing wave archive"
waves_arg=""
[ -f "$WAVES_CSV" ] && waves_arg="--waves $WAVES_CSV"
# USE_MARCONI_WAVES=1: frame tagging, the DEM wave filter and the GNSS-R wave
# allowances use the waves AT MARCONI (archive/waves_marconi.csv: camera +
# direction-converted buoy, written hourly by owg.sh) instead of the raw
# offshore buoy, which overstates sheltered south-west seas. 0 = raw buoy.
USE_MARCONI_WAVES=1
MARCONI_WAVES="$BASE/archive/waves_marconi.csv"
if [ "$USE_MARCONI_WAVES" = 1 ] && [ -f "$MARCONI_WAVES" ]; then
    waves_arg="--waves $MARCONI_WAVES"
    export CACO_LOCAL_WAVES=1
fi

# 5b. Tide gauge for GNSS-R QC. First run backfills 120 days.
gauge_days=7
[ -f "$GAUGE_CSV" ] || gauge_days=120
python3 "$BASE/fetch_tide_gauge.py" --station "$TIDE_GAUGE" --days "$gauge_days" --output "$GAUGE_CSV" >> "$LOG" 2>&1 \
    || log "WARNING: tide gauge $TIDE_GAUGE download failed; using the existing gauge archive"
qc_arg=""
[ -f "$GAUGE_CSV" ] && qc_arg="--gnssr-qc-reference $GAUGE_CSV"

# 5c. Optical wave gauge: owg.sh runs both wave models on the new c2
#     bright images (it also runs hourly from its own cron entry, and
#     emails the daily report; see owg.sh). After 5b, which has just
#     downloaded the buoy and gauge, so no second download here; its
#     lock keeps this from colliding with the hourly run.
"$BASE/owg.sh" --no-fetch --no-report >> "$LOG" 2>&1 \
    || log "WARNING: optical wave gauge (owg.sh) reported a problem; see logs/owg.log"

# 5. Match detections to measured water level, then draw the maps.
#    Both read from the ARCHIVE, not the working folders, so they see
#    the full accumulated record rather than just this run.
# Days with readings, columns 3/4/5 = YYYY MM DD: the first, and how many.
gnssr_first=""; gnssr_ndays=0
if [ -f "$GNSSR_SPLINE" ]; then
    gnssr_days=$(grep -v -e '^%' -e '^[[:space:]]*$' "$GNSSR_SPLINE" \
        | awk '{printf "%04d-%02d-%02d\n", $3, $4, $5}' | uniq)
    gnssr_first=$(echo "$gnssr_days" | head -1)
    gnssr_ndays=$(echo "$gnssr_days" | grep -c .)
fi
gnssr_known_first=""; gnssr_known_ndays=""
{ read -r gnssr_known_first gnssr_known_ndays < "$GNSSR_START_FILE"; } 2>/dev/null
if [ ! -f "$GNSSR_SPLINE" ]; then
    log "ERROR: GNSS-R file not found: $GNSSR_SPLINE"
    log "       Skipping contour extraction and maps. Check daily_gnss.sh."
elif [ -z "$gnssr_first" ]; then
    log "ERROR: GNSS-R file has no readings: $GNSSR_SPLINE"
    log "       Skipping contour extraction and maps. Check daily_gnss.sh."
elif [ -n "$gnssr_known_first" ] \
     && [[ "$gnssr_first" > "$(date -u -d "$gnssr_known_first +1 day" +%F)" ]]; then
    log "ERROR: GNSS-R now starts $gnssr_first, but the record starts $gnssr_known_first:"
    log "       rebuilding from it would drop every frame in between. Contours, maps, DEM"
    log "       and runup NOT rebuilt; the previous ones are kept. Check daily_gnss.sh on"
    log "       the GNSS side, or delete $GNSSR_START_FILE if this is intended."
elif [[ "${gnssr_known_ndays:-}" =~ ^[0-9]+$ ]] \
     && [ "$gnssr_ndays" -lt $(( gnssr_known_ndays - GNSSR_MAX_LOST_DAYS )) ]; then
    log "ERROR: GNSS-R now has readings on $gnssr_ndays days, but the record has had $gnssr_known_ndays:"
    log "       rebuilding from it would drop the frames of the missing days. Contours,"
    log "       maps, DEM and runup NOT rebuilt; the previous ones are kept. Check"
    log "       daily_gnss.sh on the GNSS side, or delete $GNSSR_START_FILE if intended."
else
    # remember the earliest start, and the most days covered (not just
    # now's, so a slow loss is still measured from the best record)
    keep_first="$gnssr_first"
    [ -n "$gnssr_known_first" ] && [[ "$gnssr_known_first" < "$gnssr_first" ]] && keep_first="$gnssr_known_first"
    keep_ndays="$gnssr_ndays"
    [[ "${gnssr_known_ndays:-}" =~ ^[0-9]+$ ]] && [ "$gnssr_known_ndays" -gt "$gnssr_ndays" ] \
        && keep_ndays="$gnssr_known_ndays"
    echo "$keep_first $keep_ndays" > "$GNSSR_START_FILE"
    # How far behind is GNSS-R? Last data row, columns 3/4/5 = YYYY MM DD.
    last_row=$(grep -v '^%' "$GNSSR_SPLINE" | tail -1)
    gnssr_last=$(echo "$last_row" | awk '{printf "%04d-%02d-%02d", $3, $4, $5}')
    if [ -n "$gnssr_last" ]; then
        lag_days=$(( ( $(date -u +%s) - $(date -u -d "$gnssr_last" +%s) ) / 86400 ))
        log "GNSS-R coverage ends $gnssr_last (${lag_days}d behind today)"
        if [ "$lag_days" -gt "$GNSSR_STALE_DAYS" ]; then
            log "WARNING: GNSS-R is ${lag_days} days behind (expected ~2). Processing may have"
            log "         stopped -- check daily_gnss.sh. Maps will keep rendering from older"
            log "         data and will NOT look broken, so this warning is the only signal."
        fi
    fi

    setup_arg=""
    [ -n "$SETUP_COEF" ] && [ -n "$waves_arg" ] && setup_arg="--setup-coef $SETUP_COEF"
    # GNSS-R QC report and plot (the contour and runup stages apply the
    # same QC themselves through $qc_arg).
    if [ -f "$GAUGE_CSV" ]; then
        python3 "$BASE/gnssr_qc.py" "$GNSSR_SPLINE" --reference "$GAUGE_CSV" \
            --output "$GNSSR_QC_CSV" --plot "$BASE/gnssr_qc_7day.png" >> "$LOG" 2>&1 \
            || log "WARNING: GNSS-R QC report failed (see above)"
    fi

    python3 "$BASE/extract_elevation_contours.py" "$GNSSR_SPLINE" \
        --processed-dir "$ARCHIVE_CSV" $waves_arg $setup_arg $qc_arg \
        --output "$CONTOURS" >> "$LOG" 2>&1
    if [ $? -ne 0 ]; then
        log "ERROR: contour extraction failed -- see above. Maps not regenerated."
    else
        # 5d. Consistency filter. MAP_CONTOURS is what the maps and
        #     georectification read: the filtered copy when the filter ran
        #     cleanly, otherwise the unfiltered file. Last run's filtered copy,
        #     report and diagnostic plots are removed first, so a failed (or
        #     disabled) run can never leave them looking current.
        MAP_CONTOURS="$CONTOURS"
        rm -f "$CONTOURS_QC" "$CONTOURS_QC.tmp" "$CONSISTENCY_REPORT" "$CONSISTENCY_REPORT.tmp" \
            "${CONSISTENCY_PLOT}_c1.png" "${CONSISTENCY_PLOT}_c2.png"
        if [ "$CONSISTENCY_ENABLE" = "1" ]; then
            qc_out=$(python3 "$BASE/waterline_consistency.py" "$CONTOURS" \
                --output "$CONTOURS_QC" --report "$CONSISTENCY_REPORT" \
                --plot "$CONSISTENCY_PLOT" --image-dir "$ARCHIVE_IMG" \
                --plot-days "$MAP_WEEK_DAYS" 2>&1)
            qc_rc=$?
            echo "$qc_out" >> "$LOG"
            if [ $qc_rc -eq 0 ] && [ -s "$CONTOURS_QC" ]; then
                MAP_CONTOURS="$CONTOURS_QC"
                log "waterline consistency: $(echo "$qc_out" | grep '^CONSISTENCY ' | tail -1 | cut -c13-)"
                # e.g. rows running AGAINST the water level in a camera's columns
                # (the detector wrong for most lines there): first line of each
                # (the full output, WARNING lines included, is in the log just
                # above; repeated here without the word, so station_status.py
                # counts each filter warning once)
                echo "$qc_out" | grep '^WARNING' | while read -r l; do log "  consistency filter: ${l#WARNING: }"; done
            else
                rm -f "$CONTOURS_QC" "$CONTOURS_QC.tmp" "$CONSISTENCY_REPORT" "$CONSISTENCY_REPORT.tmp" \
                    "${CONSISTENCY_PLOT}_c1.png" "${CONSISTENCY_PLOT}_c2.png"
                log "WARNING: waterline consistency filter failed (exit $qc_rc) -- maps, ground points"
                log "         and DEM use the UNFILTERED $(basename "$CONTOURS"). See above."
            fi
        else
            log "waterline consistency filter disabled (CONSISTENCY_ENABLE=$CONSISTENCY_ENABLE)"
        fi

        for cam in c1 c2; do
            # Rolling week -- the diagnostic view.
            python3 "$BASE/daily_elevation_map.py" "$MAP_CONTOURS" "$ARCHIVE_IMG" "$cam" \
                "$BASE/elevation_map_${cam}_${MAP_WEEK_DAYS}day.png" \
                --days "$MAP_WEEK_DAYS" >> "$LOG" 2>&1
            if [ $? -eq 0 ]; then
                log "map written: elevation_map_${cam}_${MAP_WEEK_DAYS}day.png"
            else
                log "WARNING: ${MAP_WEEK_DAYS}-day map failed for $cam (see above)"
            fi

            # Rolling month -- the change view.
            python3 "$BASE/daily_elevation_map.py" "$MAP_CONTOURS" "$ARCHIVE_IMG" "$cam" \
                "$BASE/elevation_map_${cam}_${MAP_MONTH_DAYS}day.png" \
                --days "$MAP_MONTH_DAYS" \
                --max-lines "$MAP_MONTH_MAX_LINES" >> "$LOG" 2>&1
            if [ $? -eq 0 ]; then
                log "map written: elevation_map_${cam}_${MAP_MONTH_DAYS}day.png"
            else
                log "WARNING: ${MAP_MONTH_DAYS}-day map failed for $cam (see above)"
            fi
        done

        # 6. Georectify to UTM and rebuild the DEM.
        #
        #    Runs only if the calibration is present. A missing
        #    calibration is reported rather than skipped silently,
        #    because the maps above would still be produced and the
        #    absence of a DEM would otherwise look like success.
        if [ "$DEM_ENABLE" != "1" ]; then
            log "DEM stage disabled (DEM_ENABLE=$DEM_ENABLE)"
        elif [ ! -f "$CAL/CACO05_c1_20240801_IO.yaml" ] \
          || [ ! -f "$CAL/CACO05_c2_20240801_IO.yaml" ]; then
            log "ERROR: calibration not found in $CAL -- skipping georectification and DEM."
        else
            python3 "$BASE/georectify.py" "$MAP_CONTOURS" "$GROUND" \
                --io-c1 "$CAL/CACO05_c1_20240801_IO.yaml" \
                --eo-c1 "$CAL/CACO05_c1_20251113_EO-CV.yaml" \
                --io-c2 "$CAL/CACO05_c2_20240801_IO.yaml" \
                --eo-c2 "$CAL/CACO05_c2_20251113_EO-CV.yaml" >> "$LOG" 2>&1
            if [ $? -ne 0 ]; then
                log "ERROR: georectification failed -- DEM not rebuilt."
            else
                georef_n=$(( $(wc -l < "$GROUND") - 1 ))
                log "georectified: $georef_n point(s) -> $(basename "$GROUND")"

                hs_arg=""
                [ -n "$waves_arg" ] && [ -n "$DEM_MAX_HS" ] && hs_arg="--max-hs $DEM_MAX_HS"
                [ -n "$DEM_MAX_DAY_OFFSET" ] && hs_arg="$hs_arg --max-day-offset $DEM_MAX_DAY_OFFSET"
                # dem_from_contours.py exit: 0 built, 4 no points in the window or
                # every camera-day rejected by the day test (a short page says
                # so), 3 grids written but NO page, 1 other failure. The old page
                # is removed first, so the email never attaches a stale one.
                python3 "$BASE/dem_from_contours.py" "$GROUND" "$DEM_STEM" $hs_arg \
                    --cell "$DEM_CELL" \
                    --min-points "$DEM_MIN_POINTS" \
                    --max-spread "$DEM_MAX_SPREAD" >> "$LOG" 2>&1
                dem_rc=$?
                if [ $dem_rc -eq 0 ]; then
                    log "DEM written: $(basename "$DEM_STEM")_dem.asc (+ spread, count, png)"
                elif [ $dem_rc -eq 3 ]; then
                    log "WARNING: DEM grids written but its page was NOT drawn -- no $(basename "$DEM_STEM")_dem.png this run (see above)"
                elif [ $dem_rc -eq 4 ]; then
                    log "WARNING: DEM not built: no waterline points left (none in the window, or the day test rejected every camera-day; its page says so; see above)"
                else
                    log "WARNING: DEM build failed (see above)"
                fi

                # 6b. Rolling-window DEM and week-to-week change.
                if [ -n "$DEM_WINDOW_DAYS" ]; then
                    python3 "$BASE/dem_from_contours.py" "$GROUND" "${DEM_STEM}_${DEM_WINDOW_DAYS}day" $hs_arg \
                        --last-days "$DEM_WINDOW_DAYS" --series-dir "$DEM_SERIES" \
                        --cell "$DEM_CELL" \
                        --min-points "$DEM_MIN_POINTS" \
                        --max-spread "$DEM_MAX_SPREAD" >> "$LOG" 2>&1
                    win_rc=$?
                    if [ $win_rc -eq 0 ] || [ $win_rc -eq 3 ]; then
                        if [ $win_rc -eq 0 ]; then
                            log "window DEM written: $(basename "$DEM_STEM")_${DEM_WINDOW_DAYS}day_dem.asc (+ dated copy in $DEM_SERIES)"
                        else
                            log "WARNING: window DEM grids written but its page was NOT drawn -- no $(basename "$DEM_STEM")_${DEM_WINDOW_DAYS}day_dem.png this run (see above)"
                        fi
                        # Week-to-week change. The reference week is rebuilt from the
                        # CURRENT ground file with this DEM's settings (like for like:
                        # the archived copies were built by the processing of their
                        # day); a degenerate one (too few cells, < 0.5 m of relief) is
                        # skipped for the nearest adequate one, and a change statement
                        # always says how many cells it rests on.
                        # The new DEM is the series copy THIS run wrote (named by
                        # the window's end date, from its _info.json), not whatever
                        # is newest in the folder: an earlier run's copy can be newer
                        # when the ground file's last date moves back.
                        win_end=$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1])).get("window_end") or "")' \
                            "${DEM_STEM}_${DEM_WINDOW_DAYS}day_info.json" 2>/dev/null)
                        new_dem=""
                        if [ -n "$win_end" ] && [ -f "$DEM_SERIES/dem_${win_end}_${DEM_WINDOW_DAYS}d_dem.asc" ]; then
                            new_dem="$DEM_SERIES/dem_${win_end}_${DEM_WINDOW_DAYS}d"
                        fi
                        python3 "$BASE/dem_change.py" --series "$DEM_SERIES" --days "$DEM_WINDOW_DAYS" \
                            ${new_dem:+--b "$new_dem"} --rebuild "$GROUND" -- $hs_arg \
                            --cell "$DEM_CELL" \
                            --min-points "$DEM_MIN_POINTS" \
                            --max-spread "$DEM_MAX_SPREAD" >> "$LOG" 2>&1 \
                            || log "WARNING: beach-change map failed (see above)"
                    elif [ $win_rc -eq 4 ]; then
                        log "WARNING: window DEM not built: no waterline points left in the last ${DEM_WINDOW_DAYS} days (none, or the day test rejected every camera-day; its page says so; see above)"
                    else
                        log "WARNING: window DEM build failed (see above)"
                    fi
                fi
            fi
        fi

        # 7. Runup from the archived C2 timestacks. Needs only GNSS-R
        #    and the calibration. Positions only: elevations would need
        #    an independently surveyed beach profile (the DEM is built
        #    from the waterlines, so reading setup off it is circular).
        if [ "$RUNUP_ENABLE" != "1" ]; then
            log "runup stage disabled (RUNUP_ENABLE=$RUNUP_ENABLE)"
        else
            for line in $RUNUP_LINES; do
                python3 "$BASE/runup_from_timestack.py" "$ARCHIVE_RAS" \
                    --camera c2 --line "$line" \
                    --gnssr "$GNSSR_SPLINE" --require-water-level $qc_arg \
                    --output "$BASE/runup_c2_line${line}.csv" >> "$LOG" 2>&1
                if [ $? -eq 0 ]; then
                    log "runup updated: runup_c2_line${line}.csv"
                else
                    log "WARNING: runup failed for C2 line $line (see above)"
                fi
            done
        fi
    fi
fi

# 8. Station health summary (station_status.py): also written to
#    station_status.txt; exit status 1 = warnings, 2 = alerts.
python3 "$BASE/station_status.py" --base "$BASE" --gnssr "$GNSSR_SPLINE" > "$BASE/station_status.txt" 2>&1
status_rc=$?
log "station status: $(head -1 "$BASE/station_status.txt")"
[ "$status_rc" -ge 1 ] && grep -E "^ +(WARN|ALERT)" "$BASE/station_status.txt" | while read -r l; do log "  $l"; done

log "=== run end ==="

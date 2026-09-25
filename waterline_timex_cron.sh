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
WAVE_BUOY=44008
WAVES_CSV="$BASE/archive/waves_${WAVE_BUOY}.csv"
DEM_MAX_HS=1.5

# If GNSS-R falls further behind than this, something has stopped --
# 2 days is normal, so this allows generous margin before complaining.
GNSSR_STALE_DAYS=5

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

# 5. Match detections to measured water level, then draw the maps.
#    Both read from the ARCHIVE, not the working folders, so they see
#    the full accumulated record rather than just this run.
if [ ! -f "$GNSSR_SPLINE" ]; then
    log "ERROR: GNSS-R file not found: $GNSSR_SPLINE"
    log "       Skipping contour extraction and maps. Check daily_gnss.sh."
else
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

    python3 "$BASE/extract_elevation_contours.py" "$GNSSR_SPLINE" \
        --processed-dir "$ARCHIVE_CSV" $waves_arg \
        --output "$CONTOURS" >> "$LOG" 2>&1
    if [ $? -ne 0 ]; then
        log "ERROR: contour extraction failed -- see above. Maps not regenerated."
    else
        for cam in c1 c2; do
            # Rolling week -- the diagnostic view.
            python3 "$BASE/daily_elevation_map.py" "$CONTOURS" "$ARCHIVE_IMG" "$cam" \
                "$BASE/elevation_map_${cam}_${MAP_WEEK_DAYS}day.png" \
                --days "$MAP_WEEK_DAYS" >> "$LOG" 2>&1
            if [ $? -eq 0 ]; then
                log "map written: elevation_map_${cam}_${MAP_WEEK_DAYS}day.png"
            else
                log "WARNING: ${MAP_WEEK_DAYS}-day map failed for $cam (see above)"
            fi

            # Rolling month -- the change view.
            python3 "$BASE/daily_elevation_map.py" "$CONTOURS" "$ARCHIVE_IMG" "$cam" \
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
            python3 "$BASE/georectify.py" "$CONTOURS" "$GROUND" \
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
                python3 "$BASE/dem_from_contours.py" "$GROUND" "$DEM_STEM" $hs_arg \
                    --cell "$DEM_CELL" \
                    --min-points "$DEM_MIN_POINTS" \
                    --max-spread "$DEM_MAX_SPREAD" >> "$LOG" 2>&1
                if [ $? -eq 0 ]; then
                    log "DEM written: $(basename "$DEM_STEM")_dem.asc (+ spread, count, png)"
                else
                    log "WARNING: DEM build failed (see above)"
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
                    --gnssr "$GNSSR_SPLINE" --require-water-level \
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

log "=== run end ==="

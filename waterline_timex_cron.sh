#!/bin/bash
# Timex waterline capture + archive, safe to run repeatedly
# ------------------------------------------------------------------
# Runs the detector over whatever timex images are currently present,
# then copies the results into a PERMANENT archive that survives
# cleanup.sh.
#
# WHY IT IS BUILT THIS WAY:
#   * waterline_detector_v5.py CLEARS its output folder at the start of
#     every run. So the working folder can never be the accumulator --
#     each run would destroy the previous day's detections. Instead the
#     detector writes to a scratch folder and results are copied into
#     archive/processed_timex with `cp -n` (no clobber), which
#     accumulates and is idempotent.
#   * Because it is idempotent, running it more often than necessary is
#     harmless. That matters: cleanup.sh runs at 19:45 AND 21:00 local
#     while captures continue to 20:59, so a single nightly run cannot
#     cover the day. Run it several times instead.
#   * Source images are deleted nightly, so anything not archived before
#     cleanup.sh is gone from this machine (though it may still be on
#     S3 via s3UpdateProducts_athina.sh).
#
# Add to crontab, e.g. three passes bracketing the cleanups:
#   35 19 * * * /mnt/I2Rgus_Data/waterline/waterline_timex_cron.sh
#   55 20 * * * /mnt/I2Rgus_Data/waterline/waterline_timex_cron.sh
#   30 12 * * * /mnt/I2Rgus_Data/waterline/waterline_timex_cron.sh
# (times are LOCAL, matching the rest of this crontab)

set -u

BASE=/mnt/I2Rgus_Data/waterline
SCRATCH="$BASE/processed_timex"
ARCHIVE_CSV="$BASE/archive/processed_timex"
ARCHIVE_IMG="$BASE/archive/images_timex"
LOG="$BASE/logs/timex_cron.log"

# Keep one image per day for use as a map backdrop. Set to "all" to
# archive every timex instead -- that allows full reprocessing later
# but grows by roughly 20 MB/day per camera.
BACKDROP_PATTERN="all"
GNSSR_SPLINE=/home/argus_user/GNSS/v4.1/products/refl_code/Files/usgs/usgs_spline_out.txt
CONTOURS="$BASE/contour_points_timex.csv"
MAP_DAYS=7
GNSSR_STALE_DAYS=5
mkdir -p "$SCRATCH" "$ARCHIVE_CSV" "$ARCHIVE_IMG" "$(dirname "$LOG")"

log() { echo "[$(date -u '+%Y-%m-%dT%H:%M:%SZ')] $*" >> "$LOG"; }

log "=== run start ==="

cd "$BASE" || { log "FATAL: cannot cd to $BASE"; exit 1; }

# 1. Detect. Bias correction is OFF: the corrections in the detector
#    were fitted against snap imagery and overshoot badly on timex
#    (C1 raw timex error is ~-10 px, the snap correction is -49.8 px).
python3 waterline_detector_v5.py \
    --image-suffix timex.jpg \
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

# 3. Keep a backdrop image per day before cleanup.sh removes it.
if [ "$BACKDROP_PATTERN" = "all" ]; then
    cp -pn /mnt/I2Rgus_Data/ImageProducts/*.timex.jpg "$ARCHIVE_IMG"/ 2>/dev/null
else
    cp -pn /mnt/I2Rgus_Data/ImageProducts/$BACKDROP_PATTERN "$ARCHIVE_IMG"/ 2>/dev/null
fi
img_count=$(find "$ARCHIVE_IMG" -name '*.jpg' 2>/dev/null | wc -l)
log "backdrop images in archive: $img_count"

# 4. Flag if today produced nothing -- this is the failure mode that
#    would otherwise go unnoticed, since the map still renders from
#    older archived data and looks fine.
today=$(date -u '+%Y-%m-%d')
today_epoch_prefix=$(date -u -d "$today" '+%s')
recent=$(find "$ARCHIVE_CSV" -name '*.csv' -newermt "$today" 2>/dev/null | wc -l)
if [ "$recent" -eq 0 ]; then
    log "WARNING: no detections archived with today's date ($today). If this repeats,"
    log "         the capture chain or the detector is failing and the rolling map will"
    log "         quietly go stale rather than erroring."
fi

# 5. Match detections to measured water level, then draw the rolling map.
#    Both read from the ARCHIVE, not the working folders, so they see
#    the full accumulated record rather than just this run.
if [ ! -f "$GNSSR_SPLINE" ]; then
    log "ERROR: GNSS-R file not found: $GNSSR_SPLINE"
    log "       Skipping contour extraction and map. Check daily_gnss.sh."
else
    # How far behind is GNSS-R? Last data row, column 3/4/5 = YYYY MM DD.
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
        --processed-dir "$ARCHIVE_CSV" \
        --output "$CONTOURS" >> "$LOG" 2>&1
    if [ $? -ne 0 ]; then
        log "ERROR: contour extraction failed -- see above. Map not regenerated."
    else
        for cam in c1 c2; do
            python3 "$BASE/daily_elevation_map.py" "$CONTOURS" "$ARCHIVE_IMG" "$cam" \
                "$BASE/elevation_map_${cam}_${MAP_DAYS}day.png" \
                --days "$MAP_DAYS" \
                --background-hour 12 --background-minute 0 >> "$LOG" 2>&1
            if [ $? -eq 0 ]; then
                log "map written: elevation_map_${cam}_${MAP_DAYS}day.png"
            else
                log "WARNING: map generation failed for $cam (see above)"
            fi
        done
    fi
fi

log "=== run end ==="

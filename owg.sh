#!/bin/bash
#
# owg.sh -- the optical wave gauge in one command
#
# Does everything, in order:
#   1. downloads the offshore buoy (NDBC 44008) and the Chatham tide
#      gauge (NOAA 8447435) into archive/
#   2. runs both wave models on every new c2 bright image still on disk
#      (Run C, and the sea-patch model if it is installed) and redraws
#      their 7-day plots
#   3. writes the report: text synopsis, pass/fail checks and the
#      science plots, in reports/owg/ (a dated copy of the text is kept
#      in reports/owg/daily/)
#   4. with --email, sends it -- the same way the GPS health email goes
#      (msmtp, from the station's Gmail account)
#
# Safe to run any time, as often as you like: images already measured
# are skipped, downloads merge into archives that only grow, and a lock
# stops two runs (yours and cron's) from writing at once.
#
# Usage (on the station):
#   cd /mnt/I2Rgus_Data/waterline
#   ./owg.sh                   update everything, write the report, print it
#   ./owg.sh --email           ... and email it
#   ./owg.sh --dry-run         show the email that would be sent, don't send
#   ./owg.sh --reprocess       recompute every image still on disk (after a
#                              model change), then report
#   ./owg.sh --days 14         report and plots over 14 days instead of 7
#   ./owg.sh --no-fetch        skip the downloads (offline)
#   ./owg.sh --install-cron    set up the automatic runs (below), once
#
# AUTOMATIC RUNS (--install-cron adds these to the crontab; local time):
#   40 6-20 * * *  owg.sh           hourly in daylight. Images are moved to
#                                   S3 by cleanup.sh at 19:45 and 21:00, so
#                                   an hourly run measures every frame
#                                   first; 19:40 and 20:40 are the last
#                                   chances before each cleanup.
#   20 8 * * *     owg.sh --email   the daily report, just after the GPS
#                                   health email (08:10).
#
# Log: logs/owg.log. Each run also writes logs/owg_last_run, which the
# report checks: a report saying owg.sh last ran a day ago means cron
# stopped.

set -u

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$BASE/logs"
LOG="$LOG_DIR/owg.log"
LOCK="$LOG_DIR/owg.lock"
STAMP="$LOG_DIR/owg_last_run"

# Same people as the GPS health email (station_health_mail.sh).
# Space-separated: msmtp takes each address as its own argument.
RECIPIENTS="emarsjanik@usgs.gov csherwood@usgs.gov cvolpano@contractor.usgs.gov"

WAVE_BUOY=44008
TIDE_GAUGE=8447435
WAVES_CSV="$BASE/archive/waves_${WAVE_BUOY}.csv"
GAUGE_CSV="$BASE/archive/gauge_${TIDE_GAUGE}.csv"
MODEL_C="$BASE/owg_models/owg_c2_H_current_C"
MODEL_PATCH="$BASE/owg_models/owg_c2_H_patch"
MODEL_V2="$BASE/owg_models/owg_c2_H_v2"
OUT_DIR="$BASE/reports/owg"

# Plots from other systems attached to the email (missing ones are noted,
# not fatal): the GPS wave-setup check made by the GNSS processing.
GPS_SETUP_PLOT="/home/argus_user/GNSS/v4.1/products/refl_code/Files/*/*_wave_setup_check.png"

EMAIL=false; DRY_RUN=false; FETCH=true; REPORT=true; REPROCESS=""; DAYS=7; INSTALL_CRON=false
while [ $# -gt 0 ]; do
    case "$1" in
        --email)        EMAIL=true ;;
        --dry-run)      EMAIL=true; DRY_RUN=true ;;
        --no-fetch)     FETCH=false ;;
        --no-report)    REPORT=false ;;
        --reprocess)    REPROCESS="--reprocess" ;;
        --days)         DAYS="$2"; shift ;;
        --install-cron) INSTALL_CRON=true ;;
        -h|--help)      sed -n '2,46p' "$0"; exit 0 ;;
        *) echo "Unknown option: $1 (try --help)" >&2; exit 2 ;;
    esac
    shift
done

mkdir -p "$LOG_DIR" "$BASE/archive"

if $INSTALL_CRON; then
    TAG="# owg.sh (optical wave gauge)"
    current=$(crontab -l 2>/dev/null)
    printf '%s\n' "$current" > "$LOG_DIR/crontab.before_owg_install"
    {
        printf '%s\n' "$current" | grep -vF "$TAG"
        echo "40 6-20 * * * $BASE/owg.sh $TAG: hourly update"
        echo "20 8 * * * $BASE/owg.sh --email $TAG: daily report email"
    } | sed '/^$/N;/^\n$/D' | crontab -
    echo "Crontab now (previous copy saved in $LOG_DIR/crontab.before_owg_install):"
    crontab -l | grep -F "$TAG"
    exit 0
fi

interactive=false; [ -t 1 ] && interactive=true
log() { echo "[$(date -u '+%Y-%m-%dT%H:%M:%SZ')] $*" >> "$LOG"; $interactive && echo "$*"; }

# One run at a time; wait up to 15 min for another to finish.
exec 9>"$LOCK"
if ! flock -w 900 9; then
    log "SKIPPED: another owg.sh run has held the lock for 15 minutes"
    exit 1
fi

# Run a step: output to the log, and to the screen when run by hand.
run() {
    if $interactive; then
        "$@" 2>&1 | tee -a "$LOG"
        return "${PIPESTATUS[0]}"
    fi
    "$@" >> "$LOG" 2>&1
}

log "=== owg.sh starting ==="
errors=""

if $FETCH; then
    log "--- 1. downloads: buoy $WAVE_BUOY, tide gauge $TIDE_GAUGE ---"
    run python3 "$BASE/fetch_buoy_waves.py" --station "$WAVE_BUOY" --output "$WAVES_CSV" \
        || { log "WARNING: buoy download failed; using the archive"; errors="$errors buoy"; }
    days=7; [ -f "$GAUGE_CSV" ] || days=120
    run python3 "$BASE/fetch_tide_gauge.py" --station "$TIDE_GAUGE" --days "$days" --output "$GAUGE_CSV" \
        || { log "WARNING: tide gauge download failed; using the archive"; errors="$errors gauge"; }
fi

# Camera pointing (pointing_check.py): one timex per hour per camera against
# its reference bank. Skipped for a camera until its bank exists
# (calibration/pointing/<cam>/); the report flags any movement.
for cam in c1 c2; do
    if ls "$BASE/calibration/pointing/$cam/"*.jpg >/dev/null 2>&1; then
        run python3 "$BASE/pointing_check.py" --camera "$cam" \
            || { log "WARNING: pointing check $cam failed"; errors="$errors pointing_$cam"; }
    fi
done

log "--- 2. wave models ---"
if [ -f "$MODEL_C.onnx" ]; then
    run python3 "$BASE/owg_live.py" --model "$MODEL_C" --waves-csv "$WAVES_CSV" $REPROCESS \
        --output "$BASE/archive/owg_c2_H.csv" --plot "$BASE/owg_c2_H_7day.png" --days "$DAYS" \
        --also "patch=$BASE/archive/owg_c2_H_patch.csv" \
        || { log "WARNING: Run C failed"; errors="$errors runC"; }
else
    log "Run C model not installed ($MODEL_C.onnx)"
fi
if [ -f "$MODEL_PATCH.onnx" ] && [ -f "$MODEL_PATCH.patch.json" ]; then
    run python3 "$BASE/owg_live.py" --model "$MODEL_PATCH" --waves-csv "$WAVES_CSV" $REPROCESS \
        --gauge-csv "$GAUGE_CSV" --output "$BASE/archive/owg_c2_H_patch.csv" \
        --plot "$BASE/owg_c2_H_patch_7day.png" --days "$DAYS" \
        --also "C=$BASE/archive/owg_c2_H.csv" \
        || { log "WARNING: patch model failed"; errors="$errors patch"; }
else
    log "patch model not installed ($MODEL_PATCH.onnx + .patch.json)"
fi

# Version 2 (train_owg_torch.py: pretrained ConvNeXt, colour, original frames).
# Runs alongside the others once owg_models/owg_c2_H_v2.onnx is installed;
# the report averages it with Run C into an "ensemble" reading.
if [ -f "$MODEL_V2.onnx" ]; then
    run python3 "$BASE/owg_live.py" --model "$MODEL_V2" --waves-csv "$WAVES_CSV" $REPROCESS \
        --output "$BASE/archive/owg_c2_H_v2.csv" --plot "$BASE/owg_c2_H_v2_7day.png" --days "$DAYS" \
        --also "C=$BASE/archive/owg_c2_H.csv" \
        || { log "WARNING: v2 model failed"; errors="$errors v2"; }
fi

# Best estimate of the waves off Marconi, hourly (camera models blended with
# the buoy converted to Marconi by buoy_transfer.py, once that is fitted):
# archive/waves_marconi.csv, read by the report, GNSS-IR QC and the DEM.
run python3 "$BASE/marconi_waves.py" --waves-csv "$WAVES_CSV" \
    || { log "WARNING: marconi_waves.py failed"; errors="$errors marconi_waves"; }

if [ -n "$errors" ]; then
    echo "error in:$errors ($(date -u '+%Y-%m-%d %H:%M')Z)" > "$STAMP"
else
    echo "ok ($(date -u '+%Y-%m-%d %H:%M')Z)" > "$STAMP"
fi

status=0
if $REPORT; then
    log "--- 3. report ---"
    mail_args=()
    if $EMAIL; then
        mail_args=(--email $RECIPIENTS)
        $DRY_RUN && mail_args+=(--dry-run)
    fi
    run python3 "$BASE/owg_report.py" --archive "$BASE/archive" --days "$DAYS" \
        --model-c "$MODEL_C" --model-patch "$MODEL_PATCH" --out-dir "$OUT_DIR" \
        --waves-csv "$WAVES_CSV" --gauge-csv "$GAUGE_CSV" --last-run "$STAMP" \
        --attach "$BASE/owg_c2_H_7day.png" --attach "$BASE/owg_c2_H_patch_7day.png" \
        --attach "$GPS_SETUP_PLOT" ${mail_args[@]+"${mail_args[@]}"}
    status=$?
    if [ "$status" -ge 100 ]; then
        log "ERROR: the report email could not be sent (msmtp); the report is in $OUT_DIR"
        echo "error in: email ($(date -u '+%Y-%m-%d %H:%M')Z)" > "$STAMP"
    elif [ "$status" -gt 0 ]; then
        log "report: $status check(s) failed -- see CHECKS in $OUT_DIR/owg_report.txt"
    fi
fi

log "=== owg.sh finished${errors:+ with problems in:$errors} ==="
[ -n "$errors" ] && exit 1
[ "$status" -ge 100 ] && exit 1
exit 0

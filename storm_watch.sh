#!/bin/bash
#
# storm_watch.sh -- the daily storm watch (storm_watch.py) from cron
#
# Reads the forecast for the beach once a day (Marconi Beach on CACO-05);
# when strong weather is forecast within 4 days it starts an event: a
# pre-storm DEM from the calm days before, a heads-up email, post-storm DEMs
# from the days after, and the change map (erosion / accretion) by email.
# The camera DEM maps the intertidal beach only (between the low and high
# waterlines). In alerts mode (STORM_MODE=alerts, a new station before its
# beach maps) it sends the forecast heads-up and a storm-over email only.
# Details: python3 storm_watch.py --help, docs/VALIDATION.md section 2.10,
# docs/STORM_WATCH_ALERTS_MODE.md.
#
# Usage (on the station):
#   cd /mnt/I2Rgus_Data/waterline
#   ./storm_watch.sh                  the daily run (what cron runs)
#   ./storm_watch.sh --dry-run        fetch and evaluate, print; write and send nothing
#   ./storm_watch.sh --status         the events and what each waits for
#   ./storm_watch.sh --check-settings check storm_watch_station.conf and storm_watch_settings.conf
#   ./storm_watch.sh --install-cron   add the daily line to the crontab, once
#   ./storm_watch.sh --remove-cron    take it out again (turns the storm watch off)
#   Any other option is passed to storm_watch.py (e.g. --no-email, --replay A B --out-dir D).
#   On a NUC set up by coastcam-station-setup (~/stationTools/crontab_new_station.txt
#   exists) setup owns the crontab: --install-cron and --remove-cron refuse; set
#   STORM_WATCH="yes" (or "no") in station.conf and run ./setup.sh --only 10.
#
# AUTOMATIC RUN (--install-cron adds it; the NUC's clock, America/Chicago):
#   50 7 * * *  storm_watch.sh   07:50 local = 12:50 UTC (CDT) / 13:50 UTC (CST):
#                                after owg.sh's 07:40 buoy update and the 06Z wave
#                                forecasts, before the GPS (08:10) and OWG (08:20)
#                                emails; far from the waterline cron (12:30,
#                                19:25, 20:55) and cleanup.sh (19:45, 21:00).
#
# A new station's waterline_station.conf (coastcam-station-setup writes it;
# CACO-05 has none): checked with waterline_station_check.py check-overlay and
# sourced, as the waterline cron does, so its exported settings (STATION_TZ,
# GEOTIFF_EPSG, STATION_NAME, ...) reach the storm DEM builds. A refused file
# stops the run (an ERROR line in logs/storm_watch.log, exit 1).
#
# Logs: logs/storm_watch.log (one line per run), logs/storm_watch_detail.log
# (everything the run printed), logs/storm_watch_last_run (station_status.py
# reads it). storm_watch.py takes its own lock (logs/storm_watch.lock); this
# wrapper takes none (a second flock here would conflict with it).

set -u

BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="$BASE/logs"
DETAIL="$LOG_DIR/storm_watch_detail.log"
RUNLOG="$LOG_DIR/storm_watch.log"
TAG="# storm_watch.sh (storm watch)"

mkdir -p "$LOG_DIR"

case "${1:-}" in
    -h|--help)
        sed -n '2,45p' "$0"
        exit 0 ;;
    --install-cron|--remove-cron)
        if [ -e "${HOME:-/nonexistent}/stationTools/crontab_new_station.txt" ]; then
            want=$([ "$1" = "--install-cron" ] && echo yes || echo no)
            echo "this NUC's crontab is built by coastcam-station-setup (~/stationTools/crontab_new_station.txt exists):" >&2
            echo "  not changed. Set STORM_WATCH=\"$want\" in station.conf, then run ./setup.sh --only 10" >&2
            echo "  (in the coastcam-station-setup folder): it writes the storm watch line ($TAG) itself." >&2
            exit 1
        fi
        if ! current=$(crontab -l 2>/dev/null); then
            echo "crontab -l failed: not changing the crontab (is this the station's user?)" >&2
            exit 1
        fi
        if ! printf '%s\n' "$current" | grep -q "waterline_timex_cron.sh"; then
            echo "the crontab has no waterline_timex_cron.sh line: is this the station's user? Not changed." >&2
            exit 1
        fi
        before=$(printf '%s\n' "$current" | grep -vF "$TAG" | grep -c .)
        printf '%s\n' "$current" > "$LOG_DIR/crontab.before_storm_watch_install"
        {
            printf '%s\n' "$current" | grep -vF "$TAG"
            [ "$1" = "--install-cron" ] && echo "50 7 * * * $BASE/storm_watch.sh $TAG: daily forecast check"
        } | sed '/^$/N;/^\n$/D' | crontab -
        after=$(crontab -l 2>/dev/null | grep -vF "$TAG" | grep -c .)
        if [ "$after" != "$before" ]; then
            echo "the crontab's other lines changed ($before -> $after): restore it with" >&2
            echo "  crontab $LOG_DIR/crontab.before_storm_watch_install" >&2
            exit 1
        fi
        echo "Crontab now (previous copy saved in $LOG_DIR/crontab.before_storm_watch_install):"
        crontab -l | grep -F "$TAG" || echo "  (no storm watch line: the storm watch is off)"
        exit 0 ;;
esac

# --dry-run prints the email recipients: never into a log file
dry=false
for a in "$@"; do
    case "$a" in --dry-run|--status|--climatology|--check-settings) dry=true ;; esac
done

# a new station's settings for the DEM builds (CACO-05 has no such file): checked, then sourced as the cron does
STATION_CONF="$BASE/waterline_station.conf"
if ! $dry && [ -f "$STATION_CONF" ]; then
    if ! python3 "$BASE/waterline_station_check.py" check-overlay --conf "$STATION_CONF" >> "$DETAIL" 2>&1; then
        echo "[$(date -u '+%Y-%m-%dT%H:%M:%SZ')] ERROR: waterline_station.conf refused by waterline_station_check.py check-overlay (see logs/storm_watch_detail.log): storm_watch.py not run" >> "$RUNLOG"
        echo "ERROR: waterline_station.conf refused (see logs/storm_watch_detail.log): storm_watch.py not run" >&2
        exit 1
    fi
    # shellcheck source=/dev/null
    . "$STATION_CONF"
fi

if [ -t 1 ]; then
    if $dry; then
        nice -n 10 python3 "$BASE/storm_watch.py" --base "$BASE" "$@"
        rc=$?
    else
        nice -n 10 python3 "$BASE/storm_watch.py" --base "$BASE" "$@" 2>&1 | tee -a "$DETAIL"
        rc=${PIPESTATUS[0]}
    fi
else
    if $dry; then
        nice -n 10 python3 "$BASE/storm_watch.py" --base "$BASE" "$@"
        rc=$?
    else
        nice -n 10 python3 "$BASE/storm_watch.py" --base "$BASE" "$@" >> "$DETAIL" 2>&1
        rc=$?
    fi
fi
if [ "$rc" -ne 0 ]; then
    echo "[$(date -u '+%Y-%m-%dT%H:%M:%SZ')] ERROR: storm_watch.py exited $rc (see logs/storm_watch_detail.log)" >> "$RUNLOG"
fi
exit "$rc"

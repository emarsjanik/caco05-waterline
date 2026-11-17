#!/usr/bin/env python3
"""
Tests For The Storm Watch (storm_watch.py, storm_watch.sh)
============================================================
Plain python3 (no pytest needed on the station):

    python3 test_storm_watch.py
    python3 test_storm_watch.py --real-captures /tmp/sw_raw     (on the NUC: real responses saved by
                                                                 ./storm_watch.sh --dry-run --save-raw /tmp/sw_raw)

WHAT (docs/VALIDATION.md 2.10):
  T1  settings: the cron's DEM options read from its VAR= lines (each once at column 0), the coded
      defaults without it; the DRIFT check: the cron's own window-DEM call and dem_change call, run by
      bash with a python3 shim, give exactly the options storm_watch.py builds with; the station settings
      file (unknown keys, bad values); owg.sh's RECIPIENTS counted (never printed).
  T2  HTTP and parsers on SYNTHETIC responses written from the providers' documentation
      (test_fixtures/storm_watch/): ISO durations, units, NWS grid expansion and staleness, VTEC, the alert
      filter, alert times, Open-Meteo (array, one object, nulls, a Cape Cod Bay cell, a stale run), the
      GFS-Wave bulletin (rollover, staleness, misread columns), headers (the generic User-Agent; NWS_CONTACT on
      api.weather.gov only), retries, the /points cache, the Open-Meteo 400 fallbacks.
  T3  the trigger rules (a)-(e) on synthetic forecasts, the onset, the outlook, the plausibility gate.
  T4  every source failing: exit 0, "blind", no event; the observed side still opens a storm.
  T5  the observed waves: layer order, conversions, qc, camera-only rows never hiding a storm, the end rule
      on the real 1 Feb 2025 ADCP storm.
  T6  the ground file: one pass on a camera-grouped file, the cache, the cron's lock held by another process,
      a cron run due, owg.sh's lock (no 44008 fetch without it), a truncated or half-written file (the cron's
      model-only series failing is not one), the usable-date rule.
  T7  the real ADCP winter 2024-25 (sig1000_waves_ALL.csv) as the observed waves, run day by day with a fake
      builder: the events, windows and run dates of VALIDATION 2.10's table (incl. the 25-27 Dec merge and
      the 16 Feb compound event, whose second storm opens its own event first); an event's own storm never
      listed as a wave event in its windows; again with 4 missed run days.
  T8  a whole event with the REAL dem_from_contours.py and dem_change.py on a synthetic station: planted
      erosion recovered, storm dates in no DEM, like for like after a +0.10 m reprocessing shift, no write
      outside archive/storm_watch/ and logs/, no --series-dir or --last-days.
  T9  false alarm, unverified (its pre DEM built after the wave-record wait), a missed storm, an update, a second
      storm after the first, a camera that turned (and one with no pointing log), first install, a crash while
      sending, a build that raises, the cron's DEM_ENABLE=0, msmtp retries (addresses only counted), the first
      post date, wave gaps.
  T9b a missed storm while the post-window cut before it is pending; a storm forecast while the one before is
      still rough (its own heads-up); a forecast storm merged into a compound event (not a false alarm); no
      provisional pre DEM over the storm before; a blip that does not shrink the planned window; camera-only waves
      (unverified, not a false alarm); a ground file with no rows; --status read-only and the state file restored;
      an email that cannot be composed; the give-up rules (a storm that never calms, a frozen ground file, a
      receding onset).
  T10 email: MIME types, recipients, retries, msmtp's errors without addresses; --dry-run writes nothing and
      sends nothing, also when it fetches (--save-raw then passes --real-captures).
  T11 replay: writes only under --out-dir, a hindcast, emails never sent; its refusals; the live forecast archive
      (gzipped, pruned after 730 days) read back by a replay, and a replay without a ground file.
  T12 climatology: synthetic NDBC text; the real ADCP as a perfect forecast (finding 7).
  T13 storm_watch.sh: bash -n, --install-cron / --remove-cron with a fake crontab, the ERROR line, --dry-run
      and --status output kept out of the detail log, the lock.
  T14 another station (a SYNTHETIC storm_watch_station.conf): no 44008/44013 outside the CACO-05 defaults table; the
      settings files' order (station file, then settings file: the last value wins); --check-settings (exit 1 on an
      unknown key, a value that does not parse, 'KEY = v', a missing station fact; WARN on Marconi's facts at another
      station); waterline_station.conf's DEM_MAX_HS after the cron's; no owg.sh recipients with a station file or in
      alerts mode (no email, a WARNING); STORM_BUOY / STORM_BUOY2 in the fetch, file names, GFS-Wave bulletin and
      labels; the Open-Meteo facing test; STORM_TZ through zoneinfo and through pandas (zoneinfo made unimportable).
  T15 alerts mode: a storm-over email, the event closed 72 h after the storm's end (and a forecast for its tail inside
      those hours no new event), one events.csv row, no DEM call and no DEM or NAVD88 wording; a false alarm with no
      pending pre DEM; full mode unchanged.
  T16 storm_watch.sh: --install-cron / --remove-cron refused on a setup-managed NUC (~/stationTools/
      crontab_new_station.txt, a temporary HOME); waterline_station.conf checked, sourced, its exports reach python.

HOW. Every network call goes to an injected fake; every subprocess but T8's real builds to a fake runner;
msmtp and crontab are fakes on PATH. Times are fixed (the trigger tests: 16 Nov 2026 12:50 UTC).

STATION FILES. None are read. (sig1000_waves_ALL.csv, waterline_timex_cron.sh, owg.sh and the SYNTHETIC
fixtures in test_fixtures/storm_watch/ are repository files; owg.sh's addresses are only counted.)
"""

import io
import os
import re
import csv
import sys
import json
import gzip
import time
import shutil
import tempfile
import subprocess
import contextlib
from pathlib import Path
from datetime import datetime, timezone, timedelta

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
sys.dont_write_bytecode = True
import storm_watch as sw   # noqa: E402
import dem_change          # noqa: E402

UTC = timezone.utc
FIX = HERE / "test_fixtures/storm_watch"
CRON_TEXT = (HERE / "waterline_timex_cron.sh").read_text()
NOW = datetime(2026, 11, 16, 12, 50, tzinfo=UTC)
FAILED = []


def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def LOCAL():
    return datetime(2026, 1, 1, 3, 0)      # the NUC's clock, far from every station job


def nosleep(s):
    return None


def make_base(tmp, name="base", settings='STORM_EMAIL_TO="a@example.org b@example.org"\n'):
    base = Path(tmp) / name
    (base / "archive").mkdir(parents=True)
    (base / "logs").mkdir()
    shutil.copy(str(HERE / "waterline_timex_cron.sh"), str(base))
    shutil.copy(str(HERE / "owg.sh"), str(base))
    if settings is not None:
        (base / "storm_watch_settings.conf").write_text(settings)
    (base / ".timex_cron.lock").write_text("")        # the cron's lock file, as on the station
    return base


def cfg_for(base, extra=""):
    """The settings as main() reads them: storm_watch_station.conf (only T14-T15 write one), then
    storm_watch_settings.conf (+ extra)."""
    cron, _ = sw.cron_settings(base)
    p = Path(base) / "storm_watch_settings.conf"
    if extra:
        p.write_text((p.read_text() if p.exists() else "") + extra)
    cfg, _ = sw.load_settings([Path(base) / "storm_watch_station.conf", p], cron)
    return cron, cfg


def ctx_for(base, now, mode="live", extra="", **kw):
    cron, cfg = cfg_for(base, extra)
    kw.setdefault("fetch", False)
    kw.setdefault("email", False)
    kw.setdefault("local_now", LOCAL)
    kw.setdefault("sleep", nosleep)
    return sw.Ctx(base, now, cfg, cron, mode=mode, **kw)


def ep(dt):
    return int(dt.timestamp())


def write_marconi(path, rows):
    """rows: (epoch, hs, best_from, qc, dir, tp)."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_utc", "epoch", "hs_buoy", "hs_buoy_marconi", "hs_camera", "n_frames", "hs_best", "hs_sigma",
                    "best_from", "n_models", "tp_s", "tp_camera_s", "dir_deg", "qc_flag", "qc_note"])
        for e, hs, bf, qc, d, tp in rows:
            w.writerow([datetime.fromtimestamp(e, tz=UTC).strftime("%Y-%m-%dT%H:%M:%SZ"), e, "", "", "", 0, hs, "",
                        bf, 0, tp, "", d, qc, ""])


def write_buoy(path, rows):
    """rows: (epoch, hs, dpd, mwd)."""
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_utc", "epoch", "wvht_m", "dpd_s", "apd_s", "mwd_deg"])
        for e, hs, tp, d in rows:
            w.writerow([datetime.fromtimestamp(e, tz=UTC).isoformat(), e, hs, tp, "", d])


def storm_hours(start, n_days, storm_at=None, storm_hs=(2.0, 2.6, 3.0, 3.2, 3.1, 2.8, 2.6), calm=0.9,
                tail_hours=20, tail_hs=2.0):
    """Hourly (epoch, hs) from `start` for n_days: calm, with a storm from `storm_at` (the ramp, then
    `tail_hours` at tail_hs, then calm)."""
    out = []
    for h in range(n_days * 24):
        t = start + timedelta(hours=h)
        hs = calm
        if storm_at is not None and t >= storm_at:
            k = int((t - storm_at).total_seconds() // 3600)
            if k < len(storm_hs):
                hs = storm_hs[k]
            elif k < len(storm_hs) + tail_hours:
                hs = tail_hs
        out.append((ep(t), hs))
    return out


def fake_runner_factory(calls=None, dem_rc=0):
    """A stand-in for every subprocess: dem_from_contours.py writes a 20 x 20 ramp DEM (adequate) and its
    info.json; dem_change.py a uniform -0.2 m change; dem_map / asc_to_geotiff / fetch_buoy_waves do nothing."""
    calls = calls if calls is not None else []

    def run(argv, log_path=None):
        calls.append(list(argv))
        script = Path(argv[1]).name
        if script == "dem_from_contours.py":
            stem = argv[3]
            a, b = argv[argv.index("--start-date") + 1], argv[argv.index("--end-date") + 1]
            if dem_rc:
                return dem_rc, "No georectified points matched.\n"
            Path(stem).parent.mkdir(parents=True, exist_ok=True)
            g = np.add.outer(np.linspace(-0.8, 1.0, 20), np.zeros(20))
            for suf, grid in (("_dem", g), ("_spread", np.full((20, 20), 0.1)), ("_count", np.full((20, 20), 9.0)),
                              ("_sigma", np.full((20, 20), 0.02))):
                dem_change.write_grid(f"{stem}{suf}.asc", grid, 420000.0, 4638000.0, 2.0)
            Path(f"{stem}_info.json").write_text(json.dumps({"first_date": a, "last_date": b, "days": 3,
                                                             "frames": {"c1": 8, "c2": 8},
                                                             "filters": {"camera_days_rejected": []}}))
            return 0, "NOTE: fake build\n"
        if script == "dem_change.py":
            A, B = argv[argv.index("--a") + 1], argv[argv.index("--b") + 1]
            out = Path(argv[argv.index("--output-dir") + 1])
            stem = out / ("change_" + Path(A).name.replace("dem_", "") + "_to_" + Path(B).name.replace("dem_", ""))
            d = np.full((20, 20), -0.2)
            dem_change.write_grid(f"{stem}_diff.asc", d, 420000.0, 4638000.0, 2.0)
            dem_change.write_grid(f"{stem}_sig.asc", d, 420000.0, 4638000.0, 2.0)
            return 0, "level of detection      : median 0.120 m (95%, floor 0.05 m)\n"
        if script == "fetch_buoy_waves.py":
            return 1, "offline\n"
        return 0, ""
    run.calls = calls
    return run


def hold_lock(path):
    """Another process holding a lock (util-linux flock, as the cron and owg.sh take theirs)."""
    p = subprocess.Popen(["flock", str(path), "sleep", "60"], start_new_session=True)
    time.sleep(0.5)
    return p


def drop_lock(p):
    import signal
    try:
        os.killpg(p.pid, signal.SIGKILL)
    except OSError:
        pass
    p.wait()


def tree(root, skip=()):
    """{relative path: (size, mtime)} of every file under root."""
    out = {}
    for p in Path(root).rglob("*"):
        if p.is_file() and not any(str(p).startswith(str(s)) for s in skip):
            st = p.stat()
            out[str(p.relative_to(root))] = (st.st_size, st.st_mtime_ns)
    return out


# --------------------------------------------------------------------------------------------
# T1 settings and the cron drift check
# --------------------------------------------------------------------------------------------

def shim_dir(tmp, name="bin"):
    """A python3 shim that writes its argv, NUL-separated, to $ARGV_OUT."""
    b = Path(tmp) / name
    b.mkdir(exist_ok=True)
    (b / "python3").write_text('#!/bin/bash\nprintf "%s\\0" "$@" > "$ARGV_OUT"\n')
    (b / "python3").chmod(0o755)
    return b


def cron_call_argv(tmp, block, vars_):
    """Runs a block of the cron's text with bash and the python3 shim. -> argv list."""
    b = shim_dir(tmp)
    out = Path(tmp) / "argv.bin"
    if out.exists():
        out.unlink()
    script = "".join(f'{k}="{v}"\n' for k, v in vars_.items()) + block + "\n"
    env = dict(os.environ, PATH=f"{b}:{os.environ.get('PATH', '')}", ARGV_OUT=str(out))
    subprocess.run(["bash", "-c", script], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return out.read_bytes().decode().split("\0")[:-1] if out.exists() else []


def t1():
    print("T1 settings and the cron drift check")
    cron, notes = sw.cron_settings(HERE)
    want = {"DEM_CELL": "2.0", "DEM_MIN_POINTS": "3", "DEM_MAX_SPREAD": "0.5", "DEM_MAX_HS": "1.5",
            "DEM_MAX_DAY_OFFSET": "0.15", "DEM_WINDOW_DAYS": "7", "DEM_DAY_BLOCK_LOD": "", "DEM_ENABLE": "1"}
    check(all(cron[k] == v for k, v in want.items()) and not notes,
          f"the repo cron's DEM settings read as the cron's text: {[cron[k] for k in want]}")
    with tempfile.TemporaryDirectory() as tmp:
        c2, n2 = sw.cron_settings(tmp)
        check(c2["GROUND"] == f"{tmp}/contour_points_ground.csv" and c2["LOCK"] == f"{tmp}/.timex_cron.lock"
              and c2["DEM_CELL"] == "2.0" and any("coded defaults 2.0/3/0.5/1.5/0.15/7" in n for n in n2),
              "no cron file: the coded defaults, GROUND and LOCK under the given base, and a note")
        base = make_base(tmp)
        c3, _ = sw.cron_settings(base)
        check(c3["GROUND"] == f"{base}/contour_points_ground.csv" and c3["LOCK"] == f"{base}/.timex_cron.lock"
              and c3["GNSSR_START_FILE"] == f"{base}/archive/gnssr_record_start.txt",
              "$BASE in the cron's GROUND / LOCK / GNSSR_START_FILE is the given base")
        (base / "waterline_timex_cron.sh").write_text(CRON_TEXT.replace("DEM_CELL=2.0", "DEM_CELL=two"))
        c4, n4 = sw.cron_settings(base)
        check(c4["DEM_CELL"] == "2.0" and any("DEM_CELL" in n and "not a number" in n for n in n4),
              "a non-number in the cron: the default, with a note")
    # each VAR the storm watch reads is set exactly once at column 0 (bash takes the last; we read the first)
    found = sw.var_lines(CRON_TEXT)
    once = {k: len(found.get(k, [])) for k in sw.CRON_KEYS}
    check(all(n == 1 for n in once.values()), f"each needed VAR appears once at column 0 in the cron: {once}")

    # the DRIFT check: the cron's own lines, run by bash, against dem_args()
    hs_lines = [ln for ln in CRON_TEXT.splitlines()
                if re.match(r'^\s*hs_arg=""\s*$', ln)
                or re.match(r'^\s*\[ -n "\$waves_arg" \] && \[ -n "\$DEM_MAX_HS" \] && hs_arg=', ln)
                or re.match(r'^\s*\[ -n "\$DEM_MAX_DAY_OFFSET" \] && hs_arg=', ln)]
    check(len(hs_lines) == 3, f"the cron's three hs_arg lines found ({len(hs_lines)})")
    m = re.search(r'python3 "\$BASE/dem_from_contours\.py" "\$GROUND" "\$\{DEM_STEM\}_\$\{DEM_WINDOW_DAYS\}day".*?'
                  r'>> "\$LOG" 2>&1', CRON_TEXT, re.S)
    check(m is not None, "the cron's 6b window-DEM call found")
    mc = re.search(r'python3 "\$BASE/dem_change\.py" --series.*?>> "\$LOG" 2>&1', CRON_TEXT, re.S)
    check(mc is not None, "the cron's dem_change.py call found")
    if not (m and mc and len(hs_lines) == 3):
        return
    block = "\n".join(ln.strip() for ln in hs_lines) + "\n" + m.group(0)
    cblock = "\n".join(ln.strip() for ln in hs_lines) + "\n" + mc.group(0)
    with tempfile.TemporaryDirectory() as tmp:
        bad = []
        for waves in ("--waves /x/waves_marconi.csv", ""):
            for lod in ("", "1"):
                for off in ("", "0.15"):
                    v = dict(cron)
                    v.update(BASE="/b", GROUND="/b/g.csv", DEM_STEM="/b/dem_intertidal", DEM_SERIES="/b/archive/dems",
                             LOG="/dev/null", waves_arg=waves, DEM_DAY_BLOCK_LOD=lod, DEM_MAX_DAY_OFFSET=off,
                             new_dem="")
                    argv = cron_call_argv(tmp, block, v)
                    rest = argv[3:]
                    i = rest.index("--last-days") if "--last-days" in rest else -1
                    if i >= 0 and rest[i + 1] == "7":
                        del rest[i:i + 2]
                    j = rest.index("--series-dir") if "--series-dir" in rest else -1
                    if j >= 0:
                        del rest[j:j + 2]
                    mine = sw.dem_args(v, has_hs=bool(waves))
                    if argv[:3] != ["/b/dem_from_contours.py", "/b/g.csv", "/b/dem_intertidal_7day"] or rest != mine:
                        extra = sorted(set(rest) - set(mine)) + sorted(set(mine) - set(rest))
                        bad.append((waves[:7], lod, off, extra or (rest, mine)))
                    # dem_change: options before '--', the build options after it
                    cargv = cron_call_argv(tmp, cblock, v)
                    k = cargv.index("--")
                    head, tail = cargv[1:k], cargv[k + 1:]
                    opts = sorted(x for x in head if x.startswith("--"))
                    want_opts = sorted(["--series", "--days", "--rebuild"] + (["--day-block-lod"] if lod else []))
                    if opts != want_opts or tail != [x for x in mine if x != "--day-block-se"]:
                        bad.append(("dem_change", waves[:7], lod, off, opts, tail))
        check(not bad, "the cron's window-DEM call (minus --last-days 7 --series-dir) is exactly dem_args() and its "
                       "dem_change.py call takes only --series/--days/--rebuild[/--day-block-lod] and the same build "
                       "options after '--', in all 8 cases" + (f": {bad}" if bad else ""))
    src = (HERE / "storm_watch.py").read_text()
    check('["--day-block-lod"] if ctx.cron.get("DEM_DAY_BLOCK_LOD")' in src,
          "storm_watch.py passes --day-block-lod to dem_change.py exactly when the cron does (DEM_DAY_BLOCK_LOD)")

    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "s.conf"
        p.write_text("STORM_HS_M=4.0\nSTORM_FOO=1\nSTORM_MIN_HOURS=three\nSTORM_DIR_SECTOR=340-160\n")
        cfg, notes = sw.load_settings(p, cron)
        check(cfg["STORM_HS_M"] == 4.0 and any("STORM_FOO" in n for n in notes),
              "a setting is read; an unknown key gives a note")
        check(cfg["STORM_MIN_HOURS"] == 3 and any("STORM_MIN_HOURS" in n and "default" in n for n in notes),
              "a bad value: the default, with a note")
        check(cfg["STORM_DIR_SECTOR"] == (340.0, 160.0) and bool(sw.in_sector(np.array([350.0]), cfg["STORM_DIR_SECTOR"])[0])
              and not sw.in_sector(np.array([200.0]), cfg["STORM_DIR_SECTOR"])[0], "a wrapping sector 340-160")
        cfg0, _ = sw.load_settings(None, cron)
        check(cfg0["NWS_CONTACT"] == "" and cfg0["STORM_CALM_HS"] == 1.5 and cfg0["STORM_PRE_DAYS"] == 7
              and cfg0["STORM_POST_DAYS"] == 7, "defaults: NWS_CONTACT empty; calm = DEM_MAX_HS; pre/post = DEM_WINDOW_DAYS")
    n = len(sw.owg_recipients(HERE))
    check(n == 3, f"owg.sh RECIPIENTS: {n} addresses (count only)")
    ex = sw.var_lines((HERE / "storm_watch_settings.conf.example").read_text())
    check(not ex and all(f"#{k}=" in (HERE / "storm_watch_settings.conf.example").read_text() for k in sw.SETTINGS),
          "the tracked .example sets nothing and lists every setting, commented out")
    gi = (HERE / ".gitignore").read_text().splitlines()
    check("storm_watch_settings.conf" in gi, ".gitignore keeps the station's settings file untracked")

# --------------------------------------------------------------------------------------------
# T2 HTTP and the parsers (SYNTHETIC responses written from the documentation)
# --------------------------------------------------------------------------------------------

class FakeNet:
    """An opener: routes (url substring -> list of (status, bytes) or an exception, served in turn; the last
    repeats). Records every (url, headers)."""

    def __init__(self, routes, default=None):
        self.routes = routes
        self.default = default
        self.seen = []

    def __call__(self, url, headers, timeout):
        self.seen.append((url, dict(headers)))
        for key, answers in self.routes:
            if key in url:
                a = answers.pop(0) if len(answers) > 1 else answers[0]
                if isinstance(a, Exception):
                    raise a
                return a
        if self.default is not None:
            if isinstance(self.default, Exception):
                raise self.default
            return self.default
        raise OSError("no route")


def fixture_raw():
    return {"nws_grid.json": (200, (FIX / "nws_grid_SYNTHETIC.json").read_bytes(), ""),
            "nws_alerts.json": (200, (FIX / "nws_alerts_SYNTHETIC.json").read_bytes(), ""),
            "om_marine.json": (200, (FIX / "om_marine_SYNTHETIC.json").read_bytes(), ""),
            "om_wind.json": (200, (FIX / "om_wind_SYNTHETIC.json").read_bytes(), ""),
            "gfswave_44008.bull": (200, (FIX / "gfswave_44008_SYNTHETIC.bull").read_bytes(), "")}


def alert_feature(event, ugc, vtec=None, status="Actual", mtype="Alert", onset="2026-11-17T07:00:00-05:00",
                  ends="2026-11-18T19:00:00-05:00"):
    params = {"VTEC": [vtec]} if vtec else {}
    return {"properties": {"event": event, "geocode": {"UGC": ugc} if ugc is not None else {}, "status": status,
                           "messageType": mtype, "onset": onset, "ends": ends, "expires": onset,
                           "parameters": params, "headline": event}}


def t2():
    print("T2 HTTP and the parsers (synthetic responses)")
    check([sw.parse_duration(x) for x in ("PT6H", "P1D", "P1DT12H", "P7DT21H", "PT30M")] == [6, 24, 36, 189, 1],
          "ISO durations PT6H P1D P1DT12H P7DT21H (PT30M rounds up to 1 h)")
    d, why = sw.expand_layer({"uom": "wmoUnit:km_h-1", "values": [{"validTime": "2026-11-16T12:00:00+00:00/PT2H",
                                                                   "value": 36.0}]}, "speed")
    h12 = sw.hour_of(datetime(2026, 11, 16, 12, tzinfo=UTC))
    check(not why and abs(d[h12] - 10.0) < 1e-9 and abs(d[h12 + 1] - 10.0) < 1e-9 and h12 + 2 not in d,
          "km_h-1 converted to m/s from the uom, each value filling its hours")
    d, why = sw.expand_layer({"uom": "wmoUnit:furlong_fortnight-1", "values": [{"validTime": "2026-11-16T12:00:00+00:00/PT2H",
                                                                                "value": 3}]}, "speed")
    check(not d and "not understood" in why, "an unknown uom makes the layer unusable")
    d, why = sw.expand_layer({"uom": "wmoUnit:m", "values": [{"validTime": "2026-11-16T12:00:00+00:00/PT2H",
                                                              "value": 3}]}, "speed")
    check(not d, "a length unit on a speed layer is refused")
    js = json.loads((FIX / "nws_grid_SYNTHETIC.json").read_text())
    g = sw.parse_nws_grid(js, datetime(2026, 11, 16, 12, 50, tzinfo=UTC), 96)
    w = g["wind"]
    k = int(np.flatnonzero(w["h"] == h12)[0])
    check(not g["stale"] and abs(w["ws"][k] - 18.5 / 3.6) < 1e-9 and w["wdir"][k] == 230,
          "NWS grid: wind expanded hourly, 18.5 km/h -> 5.14 m/s at 12Z, direction 230")
    wv = g["waves"]
    h19 = sw.hour_of(datetime(2026, 11, 19, 6, tzinfo=UTC))
    check(wv is not None and (h19 > wv["h"][-1] or not np.isfinite(wv["hs"][h19 - wv["h"][0]])),
          "a null value is missing (no wave height on 19 Nov)")
    g2 = sw.parse_nws_grid(js, datetime(2026, 11, 17, 12, 0, tzinfo=UTC), 96)
    check(g2["stale"], "updateTime more than 24 h old: stale")
    fc = sw.build_forecast({"nws_grid.json": (200, json.dumps(js).encode(), "")},
                           datetime(2026, 11, 17, 12, 0, tzinfo=UTC), sw.load_settings(None)[0])
    check(fc["status"]["nws_grid"].startswith("stale") and fc["wind"] is None and not fc["waves"],
          "a stale grid is shown but never used")
    m = sw.VTEC_RE.search("/O.NEW.KBOX.GL.W.0042.261117T1200Z-261119T0000Z/")
    check(m and m.group(2) == "NEW" and f"{m.group(4)}.{m.group(5)}" == "GL.W"
          and sw.vtec_time(m.group(8)) == datetime(2026, 11, 19, 0, tzinfo=UTC), "a P-VTEC string parsed")
    al = sw.parse_alerts(json.loads((FIX / "nws_alerts_SYNTHETIC.json").read_text()), ["ANZ254", "MAZ022"])
    codes = sorted(a["code"] for a in al)
    check(codes == ["CF.A", "GL.W"], f"fixture alerts: Gale Warning and Coastal Flood Watch kept; Cancel and Test dropped "
                                      f"({codes})")
    gl = [a for a in al if a["code"] == "GL.W"][0]
    check(gl["start"] == "2026-11-17T12:00:00Z" and gl["end"] == "2026-11-19T00:00:00Z",
          "alert times with -05:00 offsets read as UTC")
    feats = {"features": [
        alert_feature("Gale Warning", ["ANZ254"], "/O.CAN.KBOX.GL.W.0042.261117T1200Z-261119T0000Z/", mtype="Update"),
        alert_feature("Gale Warning", ["ANZ254"], "/O.UPG.KBOX.GL.W.0042.261117T1200Z-261119T0000Z/", mtype="Update"),
        alert_feature("High Surf Warning", ["MAZ022"], None, onset="2026-11-17T12:00:00Z", ends="2026-11-18T00:00:00Z"),
        alert_feature("Storm Warning", ["ANZ255"], "/O.NEW.KBOX.SR.W.0001.261117T1200Z-261119T0000Z/"),
        alert_feature("Coastal Flood Warning", None, "/O.NEW.KBOX.CF.W.0001.261117T1200Z-261119T0000Z/"),
        alert_feature("Gale Warning", ["ANZ254"], "/O.NEW.KBOX.GL.W.0043.261117T1200Z-261119T0000Z/", status="Exercise")]}
    al2 = sw.parse_alerts(feats, ["ANZ254", "MAZ022"])
    check(sorted(a["code"] for a in al2) == ["CF.W", "SU.W"],
          "CAN and UPG actions dropped, an Exercise dropped, another zone dropped, no-VTEC 'High Surf Warning' -> "
          "SU.W, an alert without UGC kept (the query was by zone)")
    su = [a for a in al2 if a["code"] == "SU.W"][0]
    check(su["start"] == "2026-11-17T12:00:00Z", "a 'Z' time read on Python 3.8 (fromisoformat does not take 'Z')")

    cfg = sw.load_settings(None)[0]
    pts = [cfg["STORM_MARINE_POINT"], cfg["STORM_BUOY_POINT"]]
    om = json.loads((FIX / "om_marine_SYNTHETIC.json").read_text())
    res = sw.parse_om_marine(om, pts, NOW, 96, cfg["STORM_MIN_CELL_LON"])
    check(all(isinstance(r, dict) for r in res) and res[0]["cell"] == "41.90N 69.79W",
          "Open-Meteo marine: an array of two points, each matched to its point; the grid cell recorded")
    s0 = res[0]
    k = int(np.flatnonzero(s0["h"] == sw.hour_of(datetime(2026, 11, 16, 13, tzinfo=UTC)))[0])
    check(s0["hs"][k] == s0["hs"][k - 1] and s0["tstat"] == "Tm", "3-hourly samples held 3 h; its period is a mean (Tm)")
    one = sw.parse_om_marine(om[0], pts[:1], NOW, 96, cfg["STORM_MIN_CELL_LON"])
    check(isinstance(one[0], dict), "a single object (one point) is read too")
    om_n = json.loads(json.dumps(om))
    om_n[0]["hourly"]["wave_height"][10] = None
    r = sw.parse_om_marine(om_n, pts, NOW, 96, cfg["STORM_MIN_CELL_LON"])
    check(isinstance(r[0], dict), "a null height is missing, not an error")
    om_b = json.loads(json.dumps(om))
    om_b[0]["longitude"] = -70.05
    r = sw.parse_om_marine(om_b, [(41.90, -69.90), pts[1]], NOW, 96, cfg["STORM_MIN_CELL_LON"])
    check(isinstance(r[0], str) and "Cape Cod Bay" in r[0] and isinstance(r[1], dict),
          "a Marconi cell at -70.05 (Cape Cod Bay) refused; the 44008 point kept")
    r = sw.parse_om_marine(om, pts, datetime(2026, 11, 21, 12, tzinfo=UTC), 96, cfg["STORM_MIN_CELL_LON"])
    check(all(isinstance(x, str) and "72 h" in x for x in r), "a run ending less than 72 h ahead is not fresh")
    om_x = json.loads(json.dumps(om))
    om_x[0]["hourly"]["wave_height"][12] = 99.0
    r = sw.parse_om_marine(om_x, pts, NOW, 96, cfg["STORM_MIN_CELL_LON"])
    check(isinstance(r[0], str) and "implausible" in r[0], "an Hs of 99 m makes the source unusable (plausibility gate)")
    om_u = json.loads(json.dumps(om))
    om_u[0]["hourly_units"]["wave_height"] = "ft"
    r = sw.parse_om_marine(om_u, pts, NOW, 96, cfg["STORM_MIN_CELL_LON"])
    check(isinstance(r[0], str) and "unit" in r[0], "a wave height not in m is refused")

    bull = (FIX / "gfswave_44008_SYNTHETIC.bull").read_text()
    s, cyc = sw.parse_gfs_bull(bull, NOW, 96)
    k = int(np.flatnonzero(s["h"] == sw.hour_of(datetime(2026, 11, 18, 18, tzinfo=UTC)))[0])
    check(cyc == datetime(2026, 11, 16, 6, tzinfo=UTC) and s["hs"][k] == 5.60 and s["tp"][k] == 10.5 and s["dir"][k] == 60,
          "GFS-Wave bulletin: Hst 5.60 m at 18 Nov 18Z, Tp and direction from the largest partition (10.5 s, 60)")
    roll = ("Cycle    : 20261130 18 UTC\n"
            + "".join(f"| {d:2d} {h:2d} |  {hs:.2f}  1  | * {hs:.2f}  8.0  60 |\n" for d, h, hs in
                      [(30, 18, 3.0), (30, 21, 3.4), (1, 0, 3.8), (1, 3, 3.6)] + [(1, h, 1.0) for h in range(6, 24, 3)]
                      + [(d, h, 1.0) for d in range(2, 6) for h in range(0, 24, 3)]))
    s, _ = sw.parse_gfs_bull(roll, datetime(2026, 11, 30, 20, tzinfo=UTC), 96)
    k = int(np.flatnonzero(s["h"] == sw.hour_of(datetime(2026, 12, 1, 0, tzinfo=UTC)))[0])
    check(s["hs"][k] == 3.8, "day 1 after day 30 rolls over to 1 Dec")
    try:
        sw.parse_gfs_bull(bull, datetime(2026, 11, 17, 19, tzinfo=UTC), 96)
        stale = False
    except ValueError as exc:
        stale = "stale" in str(exc)
    check(stale, "a bulletin whose cycle is more than 36 h old is stale")
    shifted = bull.replace("| 18 18 |  5.60  2  | * 4.76 10.5  60 |", "| 18 18 | 10.50  2  | * 4.76 10.5  60 |")
    try:
        sw.parse_gfs_bull(shifted, NOW, 96)
        caught = False
    except ValueError as exc:
        caught = "partitions" in str(exc)
    check(caught, "a period read into the height column (Hst 10.5 against partitions of 5.4 m) refuses the bulletin")

    # headers, retries, the /points cache, the Open-Meteo fallbacks
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, settings='NWS_CONTACT="ops@example.org"\nSTORM_EMAIL_TO="a@example.org"\n')
        pts_body = json.dumps({"properties": {"forecastGridData": "https://api.weather.gov/gridpoints/BOX/9,9",
                                              "forecast": "https://api.weather.gov/gridpoints/BOX/9,9/forecast",
                                              "forecastZone": "https://api.weather.gov/zones/forecast/MAZ022"}}).encode()
        net = FakeNet([("/points/", [(200, pts_body)]), ("/gridpoints/BOX/9,9/forecast", [(500, b"{}"), (200, b"{}")]),
                       ("/gridpoints/", [(200, (FIX / "nws_grid_SYNTHETIC.json").read_bytes())]),
                       ("/alerts/", [(200, (FIX / "nws_alerts_SYNTHETIC.json").read_bytes())]),
                       ("/zones/", [(200, b"{}")]),
                       ("marine-api", [(400, b'{"error":true,"reason":"x"}'), (400, b'{"error":true}'),
                                       (200, (FIX / "om_marine_SYNTHETIC.json").read_bytes())]),
                       ("api.open-meteo", [(200, (FIX / "om_wind_SYNTHETIC.json").read_bytes())]),
                       ("polar.ncep", [OSError("timed out")])])
        ctx = ctx_for(base, NOW, opener=net, fetch=True)
        state = sw.fresh_state(NOW)
        f = sw.fetch_forecast(ctx, state)
        uas = {u: h.get("User-Agent") for u, h in net.seen}
        check(all(ua == sw.USER_AGENT for u, ua in uas.items() if "api.weather.gov" not in u),
              "every non-NWS request: User-Agent exactly 'caco05-storm-watch (USGS coastal camera research)'")
        check(all("ops@example.org" in ua for u, ua in uas.items() if "api.weather.gov" in u),
              "NWS_CONTACT added to the User-Agent of api.weather.gov requests")
        check(not any("@" in v for u, h in net.seen if "api.weather.gov" not in u for v in h.values()),
              "no other request carries an '@' anywhere in its headers")
        check(all(h.get("Accept") == "application/geo+json" for u, h in net.seen if "api.weather.gov" in u),
              "NWS requests ask for application/geo+json")
        check(all("/points/41.89,-69.96" in u for u, _ in net.seen if "/points/" in u), "/points with <= 4 decimals")
        mar = [u for u, _ in net.seen if "marine-api" in u]
        check(len(mar) == 3 and "cell_selection=sea" in mar[0] and "forecast_days=10" in mar[1]
              and "cell_selection" not in mar[1] and "forecast_days=7" in mar[2],
              "Open-Meteo 400: retried without cell_selection (10 days), then with 7 days")
        check(len([u for u, _ in net.seen if u.endswith("/9,9/forecast")]) == 2 and f.raw["nws_forecast.json"][0] == 200,
              "a 500 is retried once and then succeeds")
        check(len([u for u, _ in net.seen if "polar.ncep" in u]) == 2 and f.raw["gfswave_44008.bull"][0] == 0,
              "no answer: retried once, then unavailable")
        check(state.get("nws_points", {}).get("forecastGridData", "").endswith("/9,9"),
              "the /points answer cached in the state (never hardcoded)")
        # a 404 from the grid: no retry of the 404, the /points cache dropped and refreshed once
        net2 = FakeNet([("/points/", [(200, pts_body)]), ("/gridpoints/BOX/9,9/forecast", [(200, b"{}")]),
                        ("/gridpoints/", [(404, b'{"title":"Not Found"}')])], default=(200, b"{}"))
        state["nws_points"]["utc"] = sw.iso_z(NOW)
        ctx2 = ctx_for(base, NOW, opener=net2, fetch=True)
        sw.fetch_forecast(ctx2, state)
        grids = [u for u, _ in net2.seen if u.endswith("/gridpoints/BOX/9,9")]
        check(len(grids) == 2 and len([u for u, _ in net2.seen if "/points/" in u]) == 1 and "nws_points" not in state,
              "a grid 404: not retried as such; /points refetched once, and the cache dropped when the grid still 404s")
        arch = Path(tmp) / "forecasts/2026-11-16"
        arch.mkdir(parents=True)
        for name, src in (("nws_grid.json", "nws_grid_SYNTHETIC.json"), ("om_marine.json", "om_marine_SYNTHETIC.json")):
            with gzip.open(str(arch / f"1250Z_{name}.gz"), "wb") as g:
                g.write((FIX / src).read_bytes())
        ctx3 = ctx_for(base, NOW, mode="replay", forecast_dir=arch, opener=FakeNet([], default=OSError("no network")))
        f3 = sw.fetch_forecast(ctx3, {})
        check(f3.raw["nws_grid.json"][0] == 200 and f3.raw["om_marine.json"][0] == 200 and not f3.calls,
              "archived responses (forecasts/<date>/<HHMMZ>_<name>.gz) are read without the network, the grid "
              "without a saved /points answer")
        st, body, why = sw.http_get("https://x/", {}, 1, opener=FakeNet([("x", [(503, b""), (503, b"")])]), sleep=nosleep)
        check(st == 503 and why.startswith("HTTP 503"), "two 5xx: unavailable with its status")
        n3 = FakeNet([("x", [(429, b""), (200, b"ok")])])
        st, body, why = sw.http_get("https://x/", {}, 1, opener=n3, sleep=nosleep)
        check(st == 200 and len(n3.seen) == 2, "a 429 is retried once")

# --------------------------------------------------------------------------------------------
# T3 the trigger
# --------------------------------------------------------------------------------------------

H_NOW = sw.hour_of(NOW)


def series(values, dirs=60.0, tps=8.0, name="om_marine offshore Marconi (test)", tstat="Tp", h0=None):
    n = len(values)
    arr = lambda v: np.full(n, float(v)) if np.isscalar(v) else np.asarray(v, float)   # noqa: E731
    return {"h": np.arange(h0 if h0 is not None else H_NOW - 12, (h0 if h0 is not None else H_NOW - 12) + n),
            "hs": arr(values), "tp": arr(tps), "dir": arr(dirs), "name": name, "tstat": tstat, "cell": "test"}


def bump(peak_time, peak, width=18.0, n=240, base_hs=1.0):
    t = np.arange(H_NOW - 12, H_NOW - 12 + n)
    p = sw.hour_of(peak_time)
    return base_hs + (peak - base_hs) * np.exp(-((t - p) / width) ** 2)


def block(start_dt, hours, hs, n=240, base_hs=1.0):
    v = np.full(n, base_hs)
    k = sw.hour_of(start_dt) - (H_NOW - 12)
    v[k:k + hours] = hs
    return v


def fc_of(waves=(), other=(), wind=None, alerts=None, group="om_marine"):
    return {"waves": list(waves), "waves_other": list(other), "wind": wind, "alerts": alerts,
            "groups": {group: list(waves)} if waves else {}, "waves_group": group if waves else None,
            "status": {}, "om_used": group == "om_marine", "notes": [], "text": []}


def alert(code, start, end, event="test alert", zones=("ANZ254",)):
    return {"code": code, "event": event, "zones": list(zones), "start": sw.iso_z(start), "end": sw.iso_z(end),
            "headline": "", "nws_headline": ""}


def t3():
    print("T3 the trigger")
    cron, _ = sw.cron_settings(HERE)
    cfg = sw.load_settings(None, cron)[0]
    pk = datetime(2026, 11, 18, 18, tzinfo=UTC)
    hs = bump(pk, 4.8)
    T = sw.evaluate_trigger(fc_of([series(hs)]), NOW, cfg)
    above = hs > 1.8
    i = int(np.argmax(above))
    on = H_NOW - 12 + i - 12
    j = len(hs) - 1 - int(np.argmax(above[::-1]))
    check(T["triggered"] and "a" in T["rules"] and T["peak"]["hs"] == 4.8 and T["peak"]["time"] == "2026-11-18T18:00:00Z",
          f"NE storm, 4.8 m at 18 Nov 18Z from 60: rule (a), peak recorded ({T['rules']})")
    check(T["onset_fc"] == sw.hour_iso(on) and T["end_fc"] == sw.hour_iso(H_NOW - 12 + j),
          f"onset = first hour > 1.8 m minus 12 h ({T['onset_fc']}); end = the last hour of that run ({T['end_fc']})")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        ctx = ctx_for(base, NOW, mode="dry")
        ctx.touched = {}
        state = sw.fresh_state(NOW)
        E = sw.handle_trigger(ctx, state, fc_of([series(hs)]), T, [], [])
        od = T["onset_fc"][:10]
        check(E is not None and E["id"] == f"{od}_fc" and E["pre"]["planned"] == [sw.dadd(od, -7), sw.dadd(od, -1)]
              and [m["kind"] for m in E["emails"]] == ["heads_up"],
              f"a new watch event {od}_fc, pre planned [onset date - 7, onset date - 1], a heads-up queued")
    T = sw.evaluate_trigger(fc_of([series(bump(pk, 4.5), dirs=220.0)]), NOW, cfg)
    check(not T["triggered"], "4.5 m from 220 (SW, sheltered): no trigger")
    T = sw.evaluate_trigger(fc_of([series(block(pk, 2, 3.6))]), NOW, cfg)
    check(not T["triggered"], "3.6 m for only 2 h: no trigger (STORM_MIN_HOURS 3)")
    T = sw.evaluate_trigger(fc_of([series(block(pk, 3, 3.6))]), NOW, cfg)
    check(T["triggered"] and T["rules"] == ["a"], "3.6 m for 3 h: rule (a)")
    T = sw.evaluate_trigger(fc_of([series(block(pk, 6, 3.2), dirs=100.0, tps=13.0)]), NOW, cfg)
    check(T["triggered"] and T["rules"] == ["d"], "3.2 m at a 13 s peak period from 100: rule (d)")
    T = sw.evaluate_trigger(fc_of([series(block(pk, 6, 3.2), dirs=100.0, tps=11.0, tstat="Tm")]), NOW, cfg)
    check(T["triggered"] and T["rules"] == ["d"], "3.2 m at an 11 s MEAN period (Open-Meteo): rule (d) at Tm >= 10 s")
    T = sw.evaluate_trigger(fc_of([series(block(pk, 6, 3.2), dirs=100.0, tps=11.0, tstat="Tp")]), NOW, cfg)
    check(not T["triggered"], "3.2 m at an 11 s PEAK period: no rule (d)")
    gl = alert("GL.W", pk - timedelta(hours=12), pk + timedelta(hours=12), "Gale Warning")
    T = sw.evaluate_trigger(fc_of([series(block(pk, 6, 3.2))], alerts=[gl]), NOW, cfg)
    check(T["triggered"] and T["rules"] == ["c"], "3.2 m with a Gale Warning: rule (c)")
    T = sw.evaluate_trigger(fc_of([series(block(pk, 6, 3.2))], alerts=[]), NOW, cfg)
    check(not T["triggered"], "3.2 m alone: no trigger")
    cfw = alert("CF.W", datetime(2026, 11, 18, 14, tzinfo=UTC), datetime(2026, 11, 18, 20, tzinfo=UTC),
                "Coastal Flood Warning", ("MAZ022",))
    T = sw.evaluate_trigger(fc_of([series(np.full(240, 1.0))], alerts=[cfw]), NOW, cfg)
    check(T["triggered"] and T["rules"] == ["b"] and T["onset_fc"] == "2026-11-18T02:00:00Z",
          "a Coastal Flood Warning alone, calm waves: rule (b), onset = the alert's start minus 12 h")
    late = alert("CF.W", NOW + timedelta(days=5), NOW + timedelta(days=6), "Coastal Flood Warning")
    T = sw.evaluate_trigger(fc_of([series(np.full(240, 1.0))], alerts=[late]), NOW, cfg)
    check(not T["triggered"], "an alert starting beyond the 96 h lead: no trigger")
    gated = dict(cfg, STORM_ALERT_WAVE_GATED=["CF.W"])
    T = sw.evaluate_trigger(fc_of([series(np.full(240, 1.0))], alerts=[cfw]), NOW, gated)
    check(not T["triggered"], "STORM_ALERT_WAVE_GATED=CF.W: a CF.W with a calm wave forecast does not trigger alone")
    T = sw.evaluate_trigger(fc_of([series(block(pk, 6, 3.2))], alerts=[cfw]), NOW, gated)
    check(T["triggered"] and T["rules"] == ["c"], "... but supports 3.2 m waves (rule (c))")
    T = sw.evaluate_trigger(fc_of([], alerts=[cfw]), NOW, gated)
    check(T["triggered"] and T["rules"] == ["b"], "... and is strong again when no wave forecast exists")
    day5 = NOW + timedelta(days=5)
    s5 = series(bump(day5, 4.8, width=6))
    f5 = fc_of([s5])
    T = sw.evaluate_trigger(f5, NOW, cfg)
    check(not T["triggered"] and any(d == day5.strftime("%Y-%m-%d") and h >= 4.7 for d, h, _ in T["outlook"]),
          "a 4.8 m peak on day 5: no trigger; in the outlook only")
    wind = {"h": np.arange(H_NOW - 12, H_NOW + 228), "ws": block(pk, 6, 18.0, base_hs=5.0),
            "gust": block(pk, 6, 24.0, base_hs=7.0), "wdir": np.full(240, 45.0), "name": "nws_grid wind"}
    T = sw.evaluate_trigger(fc_of([], wind=wind), NOW, cfg)
    check(T["triggered"] and T["rules"] == ["e"], "no wave forecast, NE wind 18 m/s for 6 h: rule (e)")
    T = sw.evaluate_trigger(fc_of([series(np.full(240, 1.0))], wind=wind), NOW, cfg)
    check(not T["triggered"], "the same wind with a (calm) wave forecast: rule (e) not used")
    other = series(block(pk, 6, 4.0), name="gfswave 44008 (test)")
    T = sw.evaluate_trigger(fc_of([series(np.full(240, 1.0))], other=[other]), NOW, cfg)
    check(not T["triggered"] and any(n.startswith("NOTE") and "gfswave" in n for n in T["notes"]),
          "the primary source calm, a secondary one at 4.0 m: a NOTE only")
    mar = series(block(pk, 6, 3.0), name="om_marine offshore Marconi (test)")
    b08 = series(block(pk, 6, 3.8), name="om_marine 44008 point (test)")
    T = sw.evaluate_trigger(fc_of([mar, b08]), NOW, cfg)
    check(T["triggered"] and all("44008 point" in b for b in T["by"]),
          "only the 44008 point reaches 3.5 m: it triggers, and says which point")
    nodir = series(block(pk, 6, 4.0), dirs=float("nan"))
    w220 = dict(wind, wdir=np.full(240, 220.0))
    T = sw.evaluate_trigger(fc_of([nodir], wind=w220), NOW, cfg)
    check(not T["triggered"], "a wave hour with no direction takes the wind's (220: out of the sector)")
    T = sw.evaluate_trigger(fc_of([nodir]), NOW, cfg)
    check(T["triggered"] and any("in the sector" in n for n in T["notes"]),
          "with neither direction: taken as in the sector, with a note")
    fc = sw.build_forecast(fixture_raw(), NOW, cfg)
    T = sw.evaluate_trigger(fc, NOW, cfg)
    check(fc["waves_group"] == "om_marine" and fc["om_used"] and T["triggered"] and "a" in T["rules"],
          f"the synthetic fixtures together: Open-Meteo first, rule (a) ({T['rules']}), onset {T['onset_fc']}")
    raw = fixture_raw()
    raw["om_marine.json"] = (500, b"", "HTTP 500")
    fc = sw.build_forecast(raw, NOW, cfg)
    check(fc["waves_group"] == "gfswave_44008" and fc["status"]["om_marine"].startswith("unavailable"),
          "Open-Meteo down: the GFS-Wave 44008 bulletin is the wave source")


# --------------------------------------------------------------------------------------------
# T4 fetch failures
# --------------------------------------------------------------------------------------------

def t4():
    print("T4 fetch failures never trigger and never crash")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        start = datetime(2026, 11, 1, tzinfo=UTC)
        rows = storm_hours(start, 20, storm_at=datetime(2026, 11, 14, 3, tzinfo=UTC))
        write_marconi(base / "archive/waves_marconi.csv", [(e, h, "blend", 1, 60, 9) for e, h in rows])
        dead = FakeNet([], default=OSError("network unreachable"))
        runner = fake_runner_factory()
        r1 = datetime(2026, 11, 12, 12, 50, tzinfo=UTC)
        ctx = ctx_for(base, r1, opener=dead, fetch=True, runner=runner)
        rc, line = sw.run_once(ctx)
        stamp = (base / "logs/storm_watch_last_run").read_text()
        st = json.loads((base / "archive/storm_watch/state.json").read_text())
        check(rc == 0 and line.startswith("blind") and stamp.startswith("blind since 2026-11-12") and not st["events"],
              f"every source down: exit 0, the run line and stamp say blind, no event ({stamp.strip()[:60]})")
        check(any(Path(a[1]).name == "fetch_buoy_waves.py" for a in runner.calls) and "44008 fetch failed" in line,
              "the stale 44008 archive was fetched (and its failure noted)")
        log = (base / "logs/storm_watch.log").read_text().splitlines()
        check(len(log) == 1 and log[0].startswith("[2026-11-12T12:50:00Z] blind"), "one line in logs/storm_watch.log")
        bad = FakeNet([], default=(200, b"not json {"))
        r2 = datetime(2026, 11, 16, 12, 50, tzinfo=UTC)
        ctx = ctx_for(base, r2, opener=bad, fetch=True, runner=runner)
        rc, line = sw.run_once(ctx)
        st = json.loads((base / "archive/storm_watch/state.json").read_text())
        ids = [e["id"] for e in st["events"]]
        check(rc == 0 and line.startswith("blind") and ids == ["2026-11-14_obs"],
              f"invalid JSON everywhere: still blind, and the observed side opens the storm of 14 Nov ({ids})")
        check(not (base / "logs/timex_cron.log").exists(), "nothing written to logs/timex_cron.log")
        raw_ok = {k: v for k, v in fixture_raw().items()}
        part = FakeNet([("marine-api", [(500, b"")]), ("api.open-meteo", [(200, raw_ok["om_wind.json"][1])]),
                        ("polar.ncep", [(200, raw_ok["gfswave_44008.bull"][1])])],
                       default=OSError("refused"))
        ctx = ctx_for(base, NOW, opener=part, fetch=True, runner=runner, extra="STORM_SOURCES=om_marine,gfswave,om_wind\n")
        rc, line = sw.run_once(ctx)
        check(rc == 0 and "om_marine unavailable" in line and "gfswave_44008 ok" in line and line.startswith("ok"),
              "a partial failure is listed in the run line; the others are used")


# --------------------------------------------------------------------------------------------
# T5 observed waves
# --------------------------------------------------------------------------------------------

def adcp_rows():
    import pandas as pd
    a = pd.read_csv(HERE / "sig1000_waves_ALL.csv", parse_dates=["time"])
    t = pd.DatetimeIndex(a["time"]).tz_localize("UTC")
    secs = ((t - pd.Timestamp("1970-01-01", tz="UTC")) // pd.Timedelta("1s")).tolist()   # any datetime resolution
    return list(zip([int(x) for x in secs], a["wh_4061"].tolist(), a["wp_peak"].tolist(), a["wvdir"].tolist()))


def t5():
    print("T5 observed waves: currency, layers, segments, the end rule")
    cron, _ = sw.cron_settings(HERE)
    cfg = sw.load_settings(None, cron)[0]
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        t0 = datetime(2026, 11, 10, tzinfo=UTC)
        e = [ep(t0 + timedelta(hours=k)) for k in range(10)]
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e[0], 1.2, "blend", 1, 60, 9), (e[1], 3.0, "buoy (raw)", 1, 60, 9),
                       (e[2], 3.0, "buoy (raw)", 1, 240, 9), (e[5], 5.0, "blend", 4, 60, 9),
                       (e[7], 1.3, "camera", 1, 60, 9)])
        write_buoy(base / "archive/waves_44013.csv", [(e[0], 3.0, 9, 60), (e[3], 3.0, 9, 60), (e[7], 2.961, 9, 60)])
        write_buoy(base / "archive/waves_44008.csv", [(e[0], 4.0, 9, 60), (e[4], 3.0, 9, 60), (e[5], 1.0, 9, 60)])
        o = sw.read_observed(base, t0 + timedelta(hours=9), cfg)
        v = lambda k: (round(float(o.prim[sw.hour_of(t0) + k - o.h0]), 3), sw.SRC_LABELS[int(o.psrc[sw.hour_of(t0) + k - o.h0])])  # noqa: E731
        check(v(0) == (1.2, "waves_marconi"), f"waves_marconi blend first, whatever the buoys say {v(0)}")
        check(v(1) == (2.49, "44008 x 0.83") and v(2) == (1.5, "44008 x 0.83"),
              f"'buoy (raw)' rows: x0.83 in the sector, x0.5 outside {v(1)} {v(2)}")
        check(v(3) == (2.632, "44013 converted"), f"no waves_marconi: 44013 converted, 0.827 Hs + 0.151 {v(3)}")
        check(v(4) == (2.49, "44008 x 0.83"), f"44008 last {v(4)}")
        check(v(5) == (0.83, "44008 x 0.83"), f"a qc_flag 4 row is dropped {v(5)}")
        k7 = sw.hour_of(t0) + 7 - o.h0
        check(abs(o.val[k7] - 2.6) < 0.01 and o.cam[k7] == 1.3,
              "44008 missing, the camera reads 1.3 m, 44013 converted 2.6 m: the calm test takes 2.6 (rough)")
        # camera-only rows never open or confirm; they never hide a storm the buoys see
        h0 = datetime(2026, 11, 1, tzinfo=UTC)
        cam = [(ep(h0 + timedelta(hours=k)), 1.3 if 100 <= k < 160 else 0.8, "camera", 1, 60, 9) for k in range(240)]
        write_marconi(base / "archive/waves_marconi.csv", cam)
        write_buoy(base / "archive/waves_44013.csv",
                   [(ep(h0 + timedelta(hours=k)), 2.961, 9, 60) for k in range(100, 160)])
        (base / "archive/waves_44008.csv").unlink()
        o = sw.read_observed(base, h0 + timedelta(hours=150), cfg)
        segs = sw.segments(o, cfg)
        check(len(segs) == 1 and segs[0]["end"] is None and segs[0]["peak"] == 2.6 and not segs[0]["camera_only"],
              "44008 stale, camera-only rows at 1.3 m, 44013 at 2.6 m: rough, no end")
        write_marconi(base / "archive/waves_marconi.csv",
                      [(ep(h0 + timedelta(hours=k)), 3.5 if 100 <= k < 110 else 0.8, "camera", 1, 60, 9)
                       for k in range(240)])
        (base / "archive/waves_44013.csv").unlink()
        o = sw.read_observed(base, h0 + timedelta(hours=200), cfg)
        segs = sw.segments(o, cfg)
        check(len(segs) == 1 and segs[0]["camera_only"] and not sw.opens(segs[0], cfg) and not sw.confirms(segs[0], cfg),
              "a camera-only rough spell (3.5 m) never opens or confirms an event")
    rows = adcp_rows()
    t1 = datetime(2025, 1, 25, tzinfo=UTC)
    sel = [(e, h) for e, h, _, _ in rows if ep(t1) <= e < ep(t1 + timedelta(days=16))]
    o = sw.Obs(sel[0][0] // 3600, len(sel))
    o.prim[:] = [h for _, h in sel]
    o.psrc[:] = 1
    s = [x for x in sw.segments(o, cfg) if sw.hour_iso(x["start"]) == "2025-02-01T11:00:00Z"]
    check(len(s) == 1 and sw.hour_iso(s[0]["end"]) == "2025-02-02T13:00:00Z" and s[0]["peak"] == 3.27,
          "the real 1 Feb 2025 storm (ADCP): onset 02-01 11Z, peak 3.27 m, end 02-02 13Z")
    gap = sw.hour_of(datetime(2025, 2, 2, 14, tzinfo=UTC)) - o.h0
    o.prim[gap:gap + 6] = np.nan
    s = [x for x in sw.segments(o, cfg) if sw.hour_iso(x["start"]) == "2025-02-01T11:00:00Z"]
    check(len(s) == 1 and sw.hour_iso(s[0]["end"]) == "2025-02-02T20:00:00Z",
          "6 hours missing (14Z-19Z) after it: a missing hour is never calm, the end moves to 02-02 20Z")


# --------------------------------------------------------------------------------------------
# T6 the ground file, the locks
# --------------------------------------------------------------------------------------------

GHEAD = ("source_file,camera,capture_time_utc,capture_epoch,pixel_column,easting_utm19,northing_utm19,"
         "tide_elevation_navd88,offshore_hs_m,offshore_tp_s,setup_correction_m,beach_elevation_navd88")


def ground_rows_simple(days, frames_per_day, hs=0.8, cams=("c1", "c2"), first_hour=11, points=3):
    """{camera: rows} of a small ground file: frames every 30 min from first_hour, `points` rows each."""
    out = {c: [] for c in cams}
    for d in days:
        for c in cams:
            n = frames_per_day(d, c) if callable(frames_per_day) else frames_per_day
            for k in range(n):
                t = datetime.fromisoformat(d + "T00:00:00+00:00") + timedelta(hours=first_hour, minutes=30 * k)
                h = hs(d, k) if callable(hs) else hs
                for p in range(points):
                    out[c].append(f"{ep(t)}.{c}.timex,{c},{t:%Y-%m-%dT%H:%M:%S}+00:00,{ep(t)},{p},{420001 + 2 * p:.2f},"
                                  f"4638001.00,0.5,{h},8.0,0.05,0.55")
    return out


def write_ground_file(path, by_cam, truncate=False):
    text = GHEAD + "\n" + "\n".join(r for c in sorted(by_cam) for r in by_cam[c]) + "\n"
    if truncate:
        text = text[:-30]
    Path(path).write_text(text)


def t6():
    print("T6 the ground file, the locks")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        days = [f"2026-11-0{d}" for d in range(1, 7)]
        frames = lambda d, c: 4 if d == "2026-11-03" else 3 if d == "2026-11-04" else 6   # noqa: E731
        hs = lambda d, k: 2.0 if (d == "2026-11-02" and k < 4) else 0.8                    # noqa: E731
        write_ground_file(base / "contour_points_ground.csv", ground_rows_simple(days, frames, hs))
        out = Path(tmp) / "w.csv"
        gi = sw.scan_ground(base / "contour_points_ground.csv", "2026-11-02", "2026-11-05", out, 1.5)
        rows = list(csv.reader(open(out)))
        check(gi["latest"] == "2026-11-06" and gi["has_hs"] and gi["cameras"] == ["c1", "c2"],
              "a camera-grouped (not time-sorted) file: latest date over the whole file, the wave column, the cameras")
        check(gi["counts"]["2026-11-02"] == {"c1": 2, "c2": 2} and gi["frames"]["2026-11-02"] == {"c1": 6, "c2": 6},
              "frames passing --max-hs counted per (date, camera)")
        check(len(rows) - 1 == 3 * 2 * (6 + 4 + 3 + 6) and all("2026-11-02" <= r[2][:10] <= "2026-11-05" for r in rows[1:]),
              "only the rows of the dates needed are copied")
        cron, cfg = cfg_for(base)
        u = sw.usable_dates(gi, "2026-11-02", "2026-11-05", cfg)
        check(u == ["2026-11-03", "2026-11-05"], f"usable: >= 4 passing frames (4 yes, 3 no, 2 of 6 rough no): {u}")
        check(sw.covered(gi, "2026-11-05") and not sw.covered(gi, "2026-11-06"), "a date is covered when the file "
                                                                                 "reaches the date after it")
        ctx = ctx_for(base, datetime(2026, 11, 7, 12, 50, tzinfo=UTC), extra="STORM_LOCK_WAIT_S=1\n")
        state = sw.fresh_state(ctx.now)
        n_scan = [0]
        real_scan = sw.scan_ground

        def counting(*a, **k):
            n_scan[0] += 1
            return real_scan(*a, **k)
        sw.scan_ground = counting
        try:
            g1 = sw.ground_info(ctx, state, "2026-11-01")
            g2 = sw.ground_info(ctx, state, "2026-11-01")
            check(g1 is not None and g2 is g1 and n_scan[0] == 1, "an unchanged ground file is not read again (mtime, size)")
            os.utime(base / "contour_points_ground.csv", (time.time() + 2, time.time() + 2))
            # the cron's lock held by another process (util-linux flock, as the cron takes it)
            p = hold_lock(base / ".timex_cron.lock")
            try:
                ctx.notes = []
                g3 = sw.ground_info(ctx, state, "2026-11-01")
                check(g3 is None and any("cron busy" in n for n in ctx.notes),
                      "the cron's lock held by another process: the ground steps skipped ('cron busy'), no crash")
            finally:
                drop_lock(p)
            g4 = sw.ground_info(ctx, state, "2026-11-01")
            check(g4 is not None and n_scan[0] == 2, "the next run reads it")
            ctx_due = ctx_for(base, ctx.now, local_now=lambda: datetime(2026, 11, 7, 12, 20))
            os.utime(base / "contour_points_ground.csv", (time.time() + 5, time.time() + 5))
            g5 = sw.ground_info(ctx_due, state, "2026-11-01")
            check(g5 is None and any("cron due" in n and "12:30" in n for n in ctx_due.notes) and n_scan[0] == 2,
                  "a cron run due at 12:30 local within 15 min: the cron's lock is not taken")
        finally:
            sw.scan_ground = real_scan
        # a truncated or half-written file
        st2 = {"ground_cache": None}
        write_ground_file(base / "contour_points_ground.csv", ground_rows_simple(days, frames, hs), truncate=True)
        ctx.notes, ctx.warnings = [], []
        g = sw.ground_info(ctx, st2, "2026-11-01")
        check(g is None and any("short" in w and w.startswith("WARNING") for w in ctx.warnings),
              "a ground file whose last line is cut: skipped with a WARNING")
        os.utime(base / "contour_points_ground.csv", (time.time() + 10, time.time() + 10))
        g = sw.ground_info(ctx, st2, "2026-11-01")
        os.utime(base / "contour_points_ground.csv", (time.time() + 20, time.time() + 20))
        g = sw.ground_info(ctx, st2, "2026-11-01")
        check(g is not None and any("accepted" in n for n in ctx.notes), "the same 3 runs in a row: accepted with a note")
        write_ground_file(base / "contour_points_ground.csv", ground_rows_simple(days, frames, hs, cams=("c1",)))
        ctx.warnings = []
        g = sw.ground_info(ctx, st2, "2026-11-01")
        check(g is None and any("c2" in w and "gone" in w for w in ctx.warnings), "a camera gone from the file: skipped")
        write_ground_file(base / "contour_points_ground.csv", ground_rows_simple(days, frames, hs))
        (base / "logs/timex_cron.log").write_text("[x] === run start ===\n[x] ERROR: georectification failed -- DEM "
                                                  "not rebuilt.\n[x] === run end ===\n")
        st3 = {"ground_cache": None}
        ctx.warnings = []
        g = sw.ground_info(ctx, st3, "2026-11-01")
        check(g is None and any("georectification failed" in w for w in ctx.warnings),
              "the cron's last run logged 'georectification failed': skipped")
        (base / "logs/timex_cron.log").write_text("[x] === run start ===\n[x] DEM written\n")
        ctx.warnings = []
        g = sw.ground_info(ctx, {"ground_cache": None}, "2026-11-01")
        check(g is None and any("did not log its end" in w for w in ctx.warnings), "a cron run that did not end: skipped")
        # the model-only DEM series (SETUP_SOURCE on) failing writes elsewhere: $GROUND was written, not suspect
        (base / "logs/timex_cron.log").write_text(
            "[x] === run start ===\n[x] georectified: 9000 point(s) -> contour_points_ground.csv\n"
            "[x] WARNING: model-only DEM series: georectification failed (see above); /b/archive/dems_model not "
            "updated\n[x] === run end ===\n")
        ctx.warnings = []
        g = sw.ground_info(ctx, {"ground_cache": None}, "2026-11-01")
        check(g is not None and not ctx.warnings,
              "the cron's model-only series logged 'georectification failed' (a WARNING, not $GROUND): the ground file "
              "is used")
        (base / "logs/timex_cron.log").write_text("[x] === run start ===\n[x] === run end ===\n")
        # owg.sh's lock held by another process
        write_marconi(base / "archive/waves_marconi.csv", [(ep(ctx.now) - 3600, 0.8, "blend", 1, 60, 9)])
        p = hold_lock(base / "logs/owg.lock")
        try:
            runner = fake_runner_factory()
            c2 = ctx_for(base, ctx.now, extra="STORM_OWG_LOCK_WAIT_S=1\n", fetch=True, runner=runner)
            o = sw.observe(c2)
            check(o.newest is not None and any("owg.lock" in n for n in c2.notes),
                  "owg.sh's lock held: the waves are read without it (noted), never waited for long")
            check(not runner.calls, "... and the stale 44008 archive is not fetched then (owg.sh writes the same file)")
        finally:
            drop_lock(p)
        runner = fake_runner_factory()
        c3 = ctx_for(base, ctx.now, fetch=True, runner=runner)
        sw.observe(c3)
        check([Path(a[1]).name for a in runner.calls] == ["fetch_buoy_waves.py"],
              "owg.sh's lock free: the stale 44008 archive is fetched under it")
        rows = ground_rows_simple(["2026-11-02", "2026-11-03"], 6, hs=lambda d, k: 2.0 if k == 0 else (0.8 if d < "2026-11-03" else ""))
        write_ground_file(Path(tmp) / "ws.csv", rows)
        ws = sw.window_stats(Path(tmp) / "ws.csv", "2026-11-02", "2026-11-03", 1.5)
        check(ws["frames"] == {"c1": 10, "c2": 10} and ws["blank_hs"] == 10 and ws["setup_m"] == 0.05
              and ws["hs_m"] == 0.8 and ws["level_range_m"] == [0.55, 0.55],
              "window statistics: frames passing --max-hs per camera, frames with no wave height (kept, counted), "
              "the median setup and Hs and the level range of the frames used")
        src = (HERE / "storm_watch.py").read_text()
        i_obs = src.index("    obs = observe(ctx)")
        i_gr = src.index("        gi = ground_info(ctx, state, need)")
        check(i_obs < i_gr and "release_lock(lockf)" in src[src.index("def observe("):src.index("def ground_info(")],
              "owg.sh's lock is released (end of observe) before the cron's is taken (ground_info): never both")

# --------------------------------------------------------------------------------------------
# T7 the real ADCP winter, day by day (fake builder)
# --------------------------------------------------------------------------------------------

def adcp_base(tmp):
    """waves_marconi.csv = the real ADCP (hs_best, qc 1) and the synthetic ground of VALIDATION 2.10: one frame per
    hour per camera at 13..20 UTC wherever the ADCP has a value, tagged with its Hs."""
    base = make_base(tmp, "adcp")
    rows = adcp_rows()
    write_marconi(base / "archive/waves_marconi.csv", [(e, round(h, 3), "blend", 1, d, round(tp, 1))
                                                         for e, h, tp, d in rows])
    counts, earliest = {}, {}
    for e, h, _, _ in rows:
        t = datetime.fromtimestamp(e, tz=UTC)
        if 13 <= t.hour <= 20 and np.isfinite(h):
            d = t.strftime("%Y-%m-%d")
            c = counts.setdefault(d, {"c1": 0, "c2": 0})
            if h <= 1.5:
                c["c1"] += 1
                c["c2"] += 1
            earliest[d] = min(earliest.get(d, 10 ** 12), e)
    return base, {"latest": max(counts), "counts": counts, "earliest": earliest, "has_hs": True}


def run_days(base, gi, d0, d1, skip=(), runner=None, lag=2):
    runner = runner or fake_runner_factory()
    when_ = {}
    d = d0
    while d <= d1:
        if d not in skip:
            now = datetime.fromisoformat(d + "T12:50:00+00:00")
            ctx = ctx_for(base, now, runner=runner, coverage_cap=sw.dadd(d, -lag))
            fc = sw.build_forecast({}, now, ctx.cfg)
            sw.run_once(ctx, ground_override=gi, fc_override=fc)
            for eid, what in ctx.touched.items():
                for w in what:
                    when_.setdefault(eid, []).append((d, w))
        d = sw.dadd(d, 1)
    st = json.loads((base / "archive/storm_watch/state.json").read_text())
    return {E["id"]: E for E in st["events"]}, when_, runner


def first_run(when_, eid, prefix):
    return next((d for d, w in when_.get(eid, []) if w.startswith(prefix)), None)


EXPECTED_ADCP = [
    # id, opened, onset, peak (m, time), end, end found, storm dates, pre (window or INADEQUATE, run),
    # post3 (window, run), post7 (window, run)
    ("2024-12-10_obs", "2024-12-12", "2024-12-10T03:00:00Z", (3.06, "2024-12-12T03:00:00Z"), "2024-12-12T08:00:00Z",
     "2024-12-13", ["2024-12-10", "2024-12-11"], ("inadequate", "2024-12-12"),
     (["2024-12-12", "2024-12-14"], "2024-12-17"), (["2024-12-12", "2024-12-18"], "2024-12-21")),
    ("2024-12-20_obs", "2024-12-21", "2024-12-20T13:00:00Z", (3.56, "2024-12-21T03:00:00Z"), "2024-12-27T13:00:00Z",
     "2024-12-28", ["2024-12-20", "2024-12-26"], (["2024-12-13", "2024-12-19"], "2024-12-22"),
     (["2024-12-27", "2024-12-29"], "2025-01-01"), (["2024-12-27", "2024-12-31"], "2025-01-03")),
    ("2025-02-01_obs", "2025-02-02", "2025-02-01T11:00:00Z", (3.27, "2025-02-01T18:00:00Z"), "2025-02-02T13:00:00Z",
     "2025-02-03", ["2025-02-01", "2025-02-01"], (["2025-01-25", "2025-01-31"], "2025-02-03"),
     (["2025-02-02", "2025-02-04"], "2025-02-07"), (["2025-02-02", "2025-02-08"], "2025-02-11")),
    ("2025-02-12_obs", "2025-02-14", "2025-02-12T16:00:00Z", (3.16, "2025-02-16T16:00:00Z"), "2025-02-17T02:00:00Z",
     "2025-02-19", ["2025-02-12", "2025-02-16"], (["2025-02-05", "2025-02-11"], "2025-02-14"),
     (["2025-02-17", "2025-02-19"], "2025-02-22"), (["2025-02-17", "2025-02-23"], "2025-02-26")),
]


def t7():
    print("T7 the real ADCP winter 2024-25, day by day (fake builder)")
    t0 = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        base, gi = adcp_base(tmp)
        ev, when_, runner = run_days(base, gi, "2024-12-10", "2025-03-10")
        check(sorted(ev) == sorted([x[0] for x in EXPECTED_ADCP] + ["2025-02-16_obs"]),
              f"4 events closed complete and the 16 Feb storm's own event, all opened from the observed waves: "
              f"{sorted(ev)}")
        for (eid, opened, onset, (pk, pkt), end, end_run, sdates, (pre, pre_run), (p3, p3_run), (p7, p7_run)) in EXPECTED_ADCP:
            E = ev.get(eid)
            if not E:
                continue
            ob = E["observed"]
            check(first_run(when_, eid, "OPENED") == opened and ob["onset"] == onset,
                  f"{eid}: opened on the {opened} run, onset {onset[:13]}Z")
            check(ob["peak_hs"] == pk and ob["peak_time"] == pkt, f"{eid}: peak {pk} m at {pkt[:13]}Z ({ob['peak_hs']} "
                                                                    f"{ob['peak_time']})")
            ends = [d for d, w in when_[eid] if w.startswith("storm over")]
            check(ob["end"] == end and ends and ends[-1] == end_run,
                  f"{eid}: end {end[:13]}Z, found on the {end_run} run ({ob['end']}, {ends})")
            check(E["storm_dates"] == sdates, f"{eid}: storm dates {sdates[0]}..{sdates[1]} ({E['storm_dates']})")
            if pre == "inadequate":
                check(E["pre"]["status"] == "inadequate" and first_run(when_, eid, "pre INADEQUATE") == pre_run,
                      f"{eid}: pre INADEQUATE on the {pre_run} run (the record starts 12-09 18Z: 0 usable dates)")
            else:
                check(E["pre"]["status"] == "final" and E["pre"]["window"] == pre and first_run(when_, eid, "pre ") == pre_run,
                      f"{eid}: pre {pre[0]}..{pre[1]} on the {pre_run} run")
            check((E["post3"] or {}).get("window") == p3 and first_run(when_, eid, "post3") == p3_run,
                  f"{eid}: post3 {p3[0]}..{p3[1]} on the {p3_run} run ({(E['post3'] or {}).get('window')})")
            check((E["post7"] or {}).get("window") == p7 and first_run(when_, eid, "post7") == p7_run
                  and E["status"] == "closed" and E["closed_reason"] == "complete",
                  f"{eid}: post7 {p7[0]}..{p7[1]} on the {p7_run} run, then closed complete")
            kinds = [m["kind"] for m in E["emails"]]
            check(kinds.count("missed_storm") == 1 and "report_post3" in kinds and "report_post7" in kinds
                  and all(m["status"] == "dry" and m["body_file"] for m in E["emails"]),
                  f"{eid}: one missed_storm email and the two reports, bodies written, nothing sent ({kinds})")
        E2 = ev.get("2024-12-20_obs", {})
        check(first_run(when_, "2024-12-20_obs", "MERGED") == "2024-12-27"
              and any("COMPOUND" in n for n in E2.get("notes", [])),
              "2024-12-20_obs: the 25-27 Dec wave event (2.29 m, 47 h) cut its post window with 2 usable dates "
              "(23, 24 Dec): MERGED on the 12-27 run")
        check(E2.get("observed", {}).get("cut_state") == "cut" and E2["observed"]["cut_date"] == "2024-12-31"
              and first_run(when_, "2024-12-20_obs", "post window CUT") == "2025-01-03",
              "... and the 1 Jan wave event (2.24 m) CUT its 7-day window at 31 Dec (4 usable dates), decided on the 01-03 run")
        c7 = (E2.get("change") or {}).get("post7") or {}
        check(any("2024-12-30" in x for x in (c7.get("waves_in_windows") or {}).get("post", [])),
              "the 30 Dec wave spell (1.70 m) inside the post7 window is listed in its change summary")
        E4 = ev.get("2025-02-12_obs", {})
        E5 = ev.get("2025-02-16_obs", {})
        check(first_run(when_, "2025-02-12_obs", "MERGED") == "2025-02-18" and E4.get("merged_from") == ["2025-02-16_obs"]
              and len(E4["observed"]["segments"]) == 2,
              "2025-02-12_obs: the 16 Feb storm (2 usable dates between) MERGED on the 02-18 run")
        check(first_run(when_, "2025-02-16_obs", "OPENED") == "2025-02-17" and E5.get("closed_reason") == "merged"
              and first_run(when_, "2025-02-16_obs", "closed: merged into 2025-02-12_obs") == "2025-02-18"
              and [m["kind"] for m in E5.get("emails", [])] == ["missed_storm", "merged"],
              "... the 16 Feb storm opened its own event on the 02-17 run, while 2025-02-12_obs's cut was still pending "
              "(its missed_storm email went out then), closed as merged on 02-18 with a short email")
        own = {"2025-02-01_obs": ("post3", ["2025-02-04T21Z-2025-02-05T06Z, peak 1.95 m"]),
               "2024-12-20_obs": ("post3", [])}
        check(all((((ev[e].get("change") or {}).get(p) or {}).get("waves_in_windows") or {}).get("post") == w
                  for e, (p, w) in own.items()),
              "an event's own storm is never listed as a wave event inside its own windows (change summary)")
        E3 = ev.get("2025-02-01_obs", {})
        check(any("2025-01-30" in f for f in E3.get("pre", {}).get("flags", []))
              and any("2025-02-09" in f for f in E4.get("pre", {}).get("flags", [])),
              "a wave event inside a pre window that cannot be avoided (30 Jan 2.09 m, 9-10 Feb 2.05 m) is flagged")
        check(not any("2025-01-20" in e for e in ev), "no event for 20 Jan 2025 (2.87 m, one hour at 2.5 m)")
        argvs = runner.calls
        check(argvs and not any("--series-dir" in a or "--last-days" in a for a in argvs),
              "no build uses --series-dir or --last-days")
    with tempfile.TemporaryDirectory() as tmp:
        base, gi = adcp_base(tmp)
        skip = ("2025-02-03", "2025-02-04", "2025-02-05", "2025-02-06")
        ev2, when2, _ = run_days(base, gi, "2025-01-26", "2025-02-14", skip=skip)
        E = ev2.get("2025-02-01_obs", {})
        check(E.get("pre", {}).get("window") == ["2025-01-25", "2025-01-31"]
              and (E.get("post3") or {}).get("window") == ["2025-02-02", "2025-02-04"]
              and (E.get("post7") or {}).get("window") == ["2025-02-02", "2025-02-08"]
              and first_run(when2, "2025-02-01_obs", "pre ") == "2025-02-07"
              and first_run(when2, "2025-02-01_obs", "post3") == "2025-02-07",
              "4 missed run days (3-6 Feb) inside the 1 Feb event: the same windows, built on the next run (7 Feb)")
    print(f"  ({time.time() - t0:.1f} s)")


# --------------------------------------------------------------------------------------------
# T8 a whole event with the REAL dem_from_contours.py and dem_change.py
# --------------------------------------------------------------------------------------------

W_TIDE = 2 * np.pi / (12.42 * 3600)
CUT_FULL, CUT_E0, CUT_E1 = 0.30, 420020.0, 420030.0


def station_rows(day, hs, eroded, storm, shift, rng):
    rows = []
    for k in range(25):                                   # 11:00-23:00 UTC every 30 min
        e_ = ep(datetime(2026, 11, day, 11, 0, tzinfo=UTC)) + 1800 * k
        tide = 0.6 + 1.2 * np.sin(W_TIDE * e_)
        t = datetime.fromtimestamp(e_, tz=UTC)
        for cam in ("c1", "c2"):
            true = tide + 0.05
            x0 = 420001.0 if cam == "c1" else 420021.0
            dz = rng.normal(0, 0.02) + (0.5 if storm else 0.0)
            for col in range(20):
                e = x0 + 2 * col
                cut = CUT_FULL * min(max((CUT_E1 - e) / (CUT_E1 - CUT_E0), 0.0), 1.0) if eroded else 0.0
                xr = int(round((true + 1.0 + cut) / 0.2))      # 0.2 m of beach per 2 m row
                rows.append((cam, f"{e_}.{cam}.timex,{cam},{t:%Y-%m-%dT%H:%M:%S}+00:00,{e_},{col},{e:.2f},"
                                  f"{4638001.0 + 2 * xr:.2f},{tide + dz + shift:.4f},{hs:.2f},8.0,0.05,"
                                  f"{true + dz + shift:.4f}"))
    return rows


def write_station_ground(path, last_day, shift=0.0):
    """Calm 1-9 Nov (Hs 0.8); a storm 10-11 Nov (Hs 3.0, lines +0.5 m); calm 12-25 Nov with the erosion ramp
    (0.30 m west of E 420020, tapering to 0 at 420030). Grouped by camera, as georectify.py writes it."""
    rng = np.random.RandomState(3)
    rows = []
    for d in range(1, last_day + 1):
        if d <= 9:
            rows += station_rows(d, 0.8, False, False, shift, rng)
        elif d <= 11:
            rows += station_rows(d, 3.0, True, True, shift, rng)
        else:
            rows += station_rows(d, 0.9, True, False, shift, rng)
    rows.sort(key=lambda r: r[0])
    Path(path).write_text(GHEAD + "\n" + "\n".join(r[1] for r in rows) + "\n")


def station_hs(h):
    """Marconi Hs at hour h since 1 Nov 00Z: onset 10 Nov 03Z, peak 3.4 m, end 12 Nov 02Z."""
    k = h - (9 * 24 + 3)
    if k < 0:
        return 0.8
    ramp = [1.6, 2.0, 2.6, 3.0, 3.2, 3.4, 3.3]
    if k < len(ramp):
        return ramp[k]
    if h < 10 * 24 + 23:
        return 3.0
    if h < 11 * 24 + 2:
        return 1.8
    return 0.9


def om_marine_json(now):
    """An Open-Meteo marine answer carrying the planned storm offshore (Marconi / 0.83), both points, hourly."""
    t0 = ep(datetime(2026, 11, 1, tzinfo=UTC))
    h0 = int(now.timestamp() // 3600) * 3600
    times, hs = [], []
    for k in range(-12, 240):
        e = h0 + 3600 * k
        times.append(datetime.fromtimestamp(e, tz=UTC).strftime("%Y-%m-%dT%H:%M"))
        hs.append(round(station_hs((e - t0) // 3600) / 0.83, 2))

    def obj(lat, lon):
        return {"latitude": lat, "longitude": lon, "hourly_units": {"wave_height": "m"},
                "hourly": {"time": times, "wave_height": hs, "wave_period": [8.0] * len(times),
                           "wave_direction": [60] * len(times)}}
    return json.dumps([obj(41.9, -69.79), obj(40.5, -69.25)]).encode()


def t8():
    print("T8 a whole event with the real dem_from_contours.py and dem_change.py (synthetic station)")
    t_start = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "station")
        (base / "archive/dems").mkdir()
        (base / "contour_points_timex.csv").write_text("sentinel\n")
        (base / "dem_intertidal_7day_dem.tif").write_bytes(b"sentinel")
        t0 = ep(datetime(2026, 11, 1, tzinfo=UTC))
        write_marconi(base / "archive/waves_marconi.csv",
                      [(t0 + 3600 * h, station_hs(h), "blend", 1, 60, 9.0) for h in range(28 * 24)])
        (base / "logs/timex_cron.log").write_text("[x] === run start ===\n[x] === run end ===\n")
        sentinels = {p: (base / p).read_bytes() for p in ("contour_points_timex.csv", "dem_intertidal_7day_dem.tif")}
        calls = []

        def recording(argv, log_path=None):
            calls.append(list(argv))
            return sw.subprocess_runner(argv, log_path)
        mine = {"contour_points_ground.csv", "archive/waves_marconi.csv"}
        before = tree(base)
        heads_up = None
        for day in range(3, 23):
            D = f"2026-11-{day:02d}"
            now = datetime.fromisoformat(D + "T12:50:00+00:00")
            write_station_ground(base / "contour_points_ground.csv", day - 2, 0.10 if day >= 18 else 0.0)
            ctx = ctx_for(base, now, runner=recording)
            raw = {"om_marine.json": (200, om_marine_json(now), "")} if day == 6 else {}
            fc = sw.build_forecast(raw, now, ctx.cfg)
            rc, line = sw.run_once(ctx, fc_override=fc)
            st = json.loads((base / "archive/storm_watch/state.json").read_text())
            if day == 6:
                E = st["events"][0] if st["events"] else {}
                check(E.get("id") == "2026-11-09_fc" and E.get("status") == "watch"
                      and (E.get("pre_provisional") or {}).get("status") == "built",
                      f"the 6 Nov forecast opens 2026-11-09_fc in WATCH with a provisional pre DEM ({line[:80]})")
                bf = evdir_file(base, E, "heads_up")
                heads_up = bf.read_text() if bf else ""
            if day == 10:
                check(st["events"][0]["status"] == "storm", "10 Nov: the observed waves confirm it (STORM)")
            if day == 13:
                check(st["events"][0]["status"] == "post", "13 Nov: the storm is over (POST)")
        pp = st["events"][0].get("pre_provisional") or {}
        check("intertidal" in heads_up and "2026-11-02 to 2026-11-08" in heads_up and "Open-Meteo.com" in heads_up
              and "pre_2026" in heads_up.split("Attachments:")[1].split("\n")[0]
              and f"{pp.get('cells')} cells" in heads_up,
              "the heads-up: intertidal note, the planned pre window, the Open-Meteo attribution (it used Open-Meteo), "
              "the provisional pre map attached")
        E = st["events"][0]
        ed = base / "archive/storm_watch/events" / E["id"]
        check(E["status"] == "closed" and E["closed_reason"] == "complete", f"closed complete ({E['history'][-1]})")
        check(E["storm_dates"] == ["2026-11-10", "2026-11-11"] and E["pre"]["window"] == ["2026-11-03", "2026-11-09"],
              "storm dates 10-11 Nov; pre final 3-9 Nov = [onset date - 7, onset date - 1]")
        info = json.loads((ed / (E["pre"]["stem"] + "_info.json")).read_text())
        check("2026-11-03" <= info["first_date"] and info["last_date"] <= "2026-11-09",
              f"the pre DEM's own dates lie inside its window ({info['first_date']}..{info['last_date']})")
        check(E["post3"]["window"] == ["2026-11-12", "2026-11-14"] and E["post7"]["window"] == ["2026-11-12", "2026-11-18"],
              "post3 12-14 Nov, post7 12-18 Nov")
        wins = []
        for a in calls:
            if Path(a[1]).name == "dem_from_contours.py":
                wins.append((a[a.index("--start-date") + 1], a[a.index("--end-date") + 1]))
        check(wins and all(not (w[0] <= "2026-11-11" and w[1] >= "2026-11-10") for w in wins),
              f"no DEM window touches the storm dates ({len(wins)} builds)")
        check(not any("--series-dir" in a or "--last-days" in a for a in calls),
              "no subprocess was given --series-dir or --last-days")
        opts = [a[a.index("--end-date") + 2:] for a in calls if Path(a[1]).name == "dem_from_contours.py"]
        want = sw.dem_args(sw.cron_settings(base)[0], True)
        check(all(o[:len(want)] == want for o in opts), f"every build carries exactly the cron's options ({' '.join(want)})")
        c3 = E["change"]["post3"]
        planted = None
        if c3.get("status") == "map":
            diff, h = dem_change.read_grid(str(ed / c3["stem"]) + "_diff.asc")
            sig, _ = dem_change.read_grid(str(ed / c3["stem"]) + "_sig.asc")
            Ec = h["xllcorner"] + h["cellsize"] * (np.arange(diff.shape[1]) + 0.5)
            full = np.broadcast_to(Ec < CUT_E0, diff.shape)
            none = np.broadcast_to(Ec >= CUT_E1, diff.shape)
            cutv = np.broadcast_to(-CUT_FULL * np.clip((CUT_E1 - Ec) / (CUT_E1 - CUT_E0), 0, 1), diff.shape)
            ok = np.isfinite(diff)
            med = float(np.median(diff[ok & full]))
            planted = float(np.sum(cutv[ok])) * h["cellsize"] ** 2
            vol = c3["erosion"]["volume_m3"]
            check(abs(med + CUT_FULL) <= 0.06, f"change post3: full-cut cells median {med:+.3f} m (planted -0.30 +/- 0.06)")
            check(abs(vol - planted) <= 0.15 * abs(planted),
                  f"erosion volume {vol:+.1f} m3 within 15 % of the planted {planted:+.1f} m3 over the compared cells")
            check(int((np.isfinite(sig) & none).sum()) == 0 and c3["accretion"]["cells"] == 0,
                  "no significant cell where nothing was planted; no accretion")
            txt = (ed / "change_post3/dem_change.txt").read_text()
            m_e = re.search(r"significant erosion\s*:\s*(\d+) cells, \d+ m2, ([+-][\d.]+) m3", txt)
            m_n = re.search(r"net \(significant cells\) : ([+-][\d.]+) m3", txt)
            check(m_e and int(m_e.group(1)) == c3["erosion"]["cells"] and abs(float(m_e.group(2)) - vol) <= 0.1
                  and m_n and abs(float(m_n.group(1)) - c3["net_m3"]) <= 0.1,
                  "change_summary.json's volumes equal dem_change.py's printed lines within 0.1 m3")
            js = json.loads((ed / "change_post3/change_summary.json").read_text())
            check(js["rebuilt"] and js["storm_dates"] == ["2026-11-10", "2026-11-11"] and js["lod_median_m"]
                  and js["pre"]["window"] == ["2026-11-03", "2026-11-09"], "change_summary.json: rebuilt pre, storm dates, LoD")
        else:
            check(False, f"change post3 made a map ({c3.get('why')})")
        c7 = E["change"]["post7"]
        if c7.get("status") == "map":
            diff, h = dem_change.read_grid(str(ed / c7["stem"]) + "_diff.asc")
            Ec = h["xllcorner"] + h["cellsize"] * (np.arange(diff.shape[1]) + 0.5)
            none = np.broadcast_to(Ec >= CUT_E1, diff.shape)
            ok = np.isfinite(diff)
            med_out = float(np.median(diff[ok & none]))
            check(abs(med_out) <= 0.03 and abs(c7["drift_rebuilt_minus_original_pre_m"] - 0.10) <= 0.02,
                  f"LIKE FOR LIKE: every row +0.10 m before post7 (a reprocessing): outside the cut the change is "
                  f"{med_out:+.3f} m; the rebuilt pre sits {c7['drift_rebuilt_minus_original_pre_m']:+.3f} m from the "
                  "original record")
        else:
            check(False, f"change post7 made a map ({c7.get('why')})")
        for part in ("pre", "post3", "post7"):
            stem = ed / E[part]["stem"]
            check(Path(str(stem) + "_map.png").exists() and Path(str(stem) + "_dem.tif").exists(),
                  f"{part}: map and GeoTIFF written in the event folder")
        check(any(p.name.endswith("_diff.tif") for p in (ed / "change_post7").iterdir())
              and any(p.suffix == ".png" for p in (ed / "change_post7").iterdir()),
              "the change map's PNG and GeoTIFFs written")
        rep = (ed / E["emails"][-1]["body_file"]).read_text()
        check("NO CHANGE MAP" not in rep and "not a sediment budget" in rep and "intertidal" in rep
              and "erosion" in rep, "the post7 report: the volumes, 'not a sediment budget', the intertidal note")
        after = tree(base)
        new = sorted(p for p in after if p not in before and p not in mine)
        changed = sorted(p for p in before if p in after and before[p] != after[p] and p not in mine)
        check(new and all(p.startswith(("archive/storm_watch/", "logs/")) for p in new) and
              all(p.startswith(("archive/storm_watch/", "logs/")) for p in changed),
              f"every file created or changed lies under archive/storm_watch/ or logs/ "
              f"({[p for p in new + changed if not p.startswith(('archive/storm_watch/', 'logs/'))][:5]})")
        check(all((base / p).read_bytes() == b for p, b in sentinels.items()) and not any((base / "archive/dems").iterdir())
              and not (base / "archive/dems_model").exists(),
              "contour_points_timex.csv and dem_intertidal_7day_dem.tif unchanged; archive/dems empty; no archive/dems_model")
        check(not (base / "archive/storm_watch/work").exists(), "the working ground copy is removed once no event is open")
    print(f"  ({time.time() - t_start:.1f} s)")


def evdir_file(base, E, kind):
    for m in E.get("emails", []):
        if m["kind"] == kind and m.get("body_file"):
            return base / "archive/storm_watch/events" / E["id"] / m["body_file"]
    return None

# --------------------------------------------------------------------------------------------
# T9 false alarm, unverified, a missed storm, an update, a second storm, first install, a crash while sending
# --------------------------------------------------------------------------------------------

def simple_ground(days):
    counts = {d: {"c1": 8, "c2": 8} for d in days}
    earliest = {d: ep(datetime.fromisoformat(d + "T13:00:00+00:00")) for d in days}
    return {"latest": days[-1], "counts": counts, "earliest": earliest, "has_hs": True}


def storm_fc(ctx, peak_time, peak=4.8, alerts=None):
    """A forecast with an om_marine storm peaking at peak_time (relative to ctx.now)."""
    h0 = ctx.now_h - 12
    t = np.arange(h0, h0 + 240)
    hs = 1.0 + (peak - 1.0) * np.exp(-((t - sw.hour_of(peak_time)) / 18.0) ** 2)
    s = {"h": t, "hs": hs, "tp": np.full(240, 9.0), "dir": np.full(240, 60.0), "name": "om_marine offshore Marconi (test)",
         "tstat": "Tm", "cell": "test"}
    return fc_of([s], alerts=alerts)


def run_day(base, d, gi, fc_fn=None, runner=None, lag=2, mailer=None, email=False, **kw):
    now = datetime.fromisoformat(d + "T12:50:00+00:00")
    ctx = ctx_for(base, now, runner=runner or fake_runner_factory(), coverage_cap=sw.dadd(d, -lag), mailer=mailer,
                  email=email, **kw)
    fc = fc_fn(ctx) if fc_fn else sw.build_forecast({}, now, ctx.cfg)
    rc, line = sw.run_once(ctx, ground_override=gi, fc_override=fc)
    st = json.loads((base / "archive/storm_watch/state.json").read_text())
    return {E["id"]: E for E in st["events"]}, ctx, line


def t9():
    print("T9 false alarm, unverified, missed storm, update, a second storm, first install, a crash while sending")
    days = [f"2026-11-{d:02d}" for d in range(1, 31)]
    gi = simple_ground(days)
    start = datetime(2026, 10, 20, tzinfo=UTC)
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "fa")
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(start, 45, calm=1.0)])
        pk = datetime(2026, 11, 13, 18, tzinfo=UTC)
        ev, _, _ = run_day(base, "2026-11-10", gi, lambda c: storm_fc(c, pk))
        E = list(ev.values())[0]
        check(E["status"] == "watch" and [m["kind"] for m in E["emails"]] == ["heads_up"],
              "a forecast storm opens a WATCH event with one heads-up")
        ev, _, _ = run_day(base, "2026-11-11", gi, lambda c: storm_fc(c, pk + timedelta(hours=6)))
        E = list(ev.values())[0]
        check(len(ev) == 1 and [m["kind"] for m in E["emails"]] == ["heads_up"] and len(E["forecast"]["runs"]) == 2,
              "the next run's trigger for the same storm updates the event: no second heads-up")
        ev, _, _ = run_day(base, "2026-11-12", gi, lambda c: storm_fc(c, pk, peak=6.0))
        ev, _, _ = run_day(base, "2026-11-12", gi, lambda c: storm_fc(c, pk, peak=7.2))
        E = list(ev.values())[0]
        check([m["kind"] for m in E["emails"]].count("update") == 1,
              "the peak up by >= 1 m: an update email, at most one a day")
        end_fc = sw.parse_iso(E["forecast"]["end_fc"])
        d_close = None
        for d in [f"2026-11-{x:02d}" for x in range(13, 24)]:
            ev, ctx, _ = run_day(base, d, gi)
            E = list(ev.values())[0]
            if E["status"] == "closed" and d_close is None:
                d_close = d
        check(E["closed_reason"] == "false_alarm" and d_close is not None
              and datetime.fromisoformat(d_close + "T12:50:00+00:00") > end_fc + timedelta(hours=48)
              and [m["kind"] for m in E["emails"]].count("false_alarm") == 1,
              f"waves < 2.5 m and no trigger: closed FALSE ALARM after end_fc + 48 h ({d_close}), one email")
        check(E["pre"]["status"] == "final" and E["pre"]["window"] == E["pre"]["planned"] and not E["pending"],
              f"the false alarm's pre DEM built from its planned window once covered ({E['pre'].get('window')})")
        body = (base / "archive/storm_watch/events" / E["id"] / [m for m in E["emails"] if m["kind"] == "false_alarm"][0]
                ["body_file"]).read_text()
        check("FALSE ALARM" in body and "Open-Meteo.com" in body and "intertidal" in body,
              "the false-alarm email: one paragraph, the attribution (Open-Meteo was used), the intertidal note")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "unv")
        rows = storm_hours(start, 22, calm=1.0)
        write_marconi(base / "archive/waves_marconi.csv", [(e, h, "blend", 1, 60, 9) for e, h in rows])
        pk = datetime(2026, 11, 13, 18, tzinfo=UTC)
        run_day(base, "2026-11-10", gi, lambda c: storm_fc(c, pk))
        seen = {}
        for d in [f"2026-11-{x:02d}" for x in range(11, 27)]:
            ev, _, _ = run_day(base, d, gi)
            E = list(ev.values())[0]
            seen[d] = (E["status"], E["pre"].get("status"), len(E["pending"]))
        check(E["closed_reason"] == "unverified" and "unverified" in [m["kind"] for m in E["emails"]],
              "the wave record stops before the storm (< 50 % of its hours): closed UNVERIFIED after end_fc + 7 d")
        check(seen["2026-11-22"] == ("closed", "planned", 1) and seen["2026-11-24"] == ("closed", "planned", 1)
              and seen["2026-11-25"] == ("closed", "final", 0) and E["pre"]["window"] == ["2026-11-05", "2026-11-11"]
              and any("gaps on 2026-11-11" in f for f in E["pre"].get("flags") or []),
              f"... its pre DEM (the record has no wave on 11 Nov) stays pending while the wave record is waited for, "
              f"and is built with a flag 3 days later ({seen['2026-11-22']}, {seen['2026-11-25']})")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "missed")
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(start, 45, storm_at=storm)])
        run_day(base, "2026-11-12", gi)
        ev, _, _ = run_day(base, "2026-11-15", gi)
        E = ev.get("2026-11-14_obs")
        check(E is not None and E["status"] in ("storm", "post") and E["forecast"] is None
              and [m["kind"] for m in E["emails"]] == ["missed_storm"],
              f"no forecast, a 3.2 m storm: opened from the observations as 2026-11-14_obs ({list(ev)})")
        for d in ("2026-11-16", "2026-11-17"):
            ev, _, _ = run_day(base, d, gi)
        E = ev["2026-11-14_obs"]
        check(E["pre"]["status"] == "final" and E["pre"]["window"] == ["2026-11-07", "2026-11-13"],
              "its pre DEM built after the fact from the saved waterline points (7-13 Nov)")
        # a camera that turned in the post window (pointing_daily.py's log): each camera's change on its own
        log = {"days": {"2026-11-16": {"status": "MOVED", "d_deg": [-1.2, 0.0, 0.0]}}}
        (base / "archive/pointing_daily_c1.json").write_text(json.dumps(log))     # no pointing_daily_c2.json
        for d in ("2026-11-18", "2026-11-19", "2026-11-20"):
            ev, _, _ = run_day(base, d, gi)
        c3 = ev["2026-11-14_obs"]["change"].get("post3") or {}
        rep = evdir_file(base, ev["2026-11-14_obs"], "report_post3")
        check((c3.get("per_camera") or {}).get("turned") == ["c1"] and (c3["per_camera"].get("c1") or {}).get("cells_compared")
              and rep is not None and "CAMERA POINTING: c1 turned on 2026-11-16 MOVED" in rep.read_text()
              and "pointing of c1 not yet logged for" in rep.read_text(),
              "c1 MOVED on a post date (not in the pre window): the change of each camera on its own, flagged in the "
              "report; dates the pointing log has not reached are named")
        check("no pointing log for c2 (archive/pointing_daily_c2.json)" in rep.read_text()
              and "pointing of c2 not yet logged" not in rep.read_text(),
              "no pointing log at all for c2: the report says so (not 'not yet logged')")
        # a second storm forecast after the first one's end: its own watch event, with a heads-up
        pk2 = datetime(2026, 11, 19, 18, tzinfo=UTC)
        ev, _, _ = run_day(base, "2026-11-17", gi, lambda c: storm_fc(c, pk2))
        new = [e for e in ev.values() if e["kind"] == "fc"]
        check(len(new) == 1 and new[0]["prev_event"] == "2026-11-14_obs" and new[0]["emails"][0]["kind"] == "heads_up",
              "a storm forecast for 2 days after the first one's end: a new WATCH event (prev_event set), with a heads-up")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "install")
        storm = datetime(2026, 11, 5, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(start, 45, storm_at=storm)])
        ev, _, _ = run_day(base, "2026-11-12", gi)
        check(not ev, "first install: a storm that ended before it opens nothing and sends nothing")
        st = json.loads((base / "archive/storm_watch/state.json").read_text())
        check(len(st["seen_segments"]) == 1, "... it is remembered as seen")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "crash")
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(start, 45, storm_at=storm)])
        run_day(base, "2026-11-12", gi)
        sent = []

        def crash_after_send(subj, to, body, att):
            sent.append(subj)
            raise KeyboardInterrupt("killed after msmtp")
        try:
            run_day(base, "2026-11-15", gi, mailer=crash_after_send, email=True)
        except KeyboardInterrupt:
            pass
        st = json.loads((base / "archive/storm_watch/state.json").read_text())
        m = st["events"][0]["emails"][0]
        check(len(sent) == 1 and m["status"] == "sending", "the email is saved as 'sending' before msmtp runs")
        ev, ctx, line = run_day(base, "2026-11-16", gi, mailer=lambda *a: sent.append("again") or 0, email=True)
        E = ev["2026-11-14_obs"]
        check(len(sent) == 1 and E["emails"][0]["status"] == "unknown" and "not resent" in line,
              "a run stopped while sending: the next run marks it 'unknown' and never sends it twice")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "broken")
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(start, 45, storm_at=storm)])
        run_day(base, "2026-11-12", gi)

        def broken(argv, log_path=None):
            raise RuntimeError("a build that broke")
        with contextlib.redirect_stdout(io.StringIO()):
            ev, _, line = run_day(base, "2026-11-17", gi, runner=broken)
        stamp = (base / "logs/storm_watch_last_run").read_text()
        check("ERROR: event 2026-11-14_obs: RuntimeError" in line and stamp.startswith("error in: 2026-11-14_obs")
              and ev["2026-11-14_obs"]["pre"]["status"] == "planned",
              "a build that raises: the run goes on (exit 0), its line says ERROR with the event, the stamp 'error in'")
        ev, _, line = run_day(base, "2026-11-18", gi)
        check(ev["2026-11-14_obs"]["pre"]["status"] == "final" and "ERROR" not in line
              and (base / "logs/storm_watch_last_run").read_text().startswith(("ok", "blind")),
              "... and the next run builds it")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "nodem")
        (base / "waterline_timex_cron.sh").write_text(CRON_TEXT.replace("\nDEM_ENABLE=1\n", "\nDEM_ENABLE=0\n"))
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(start, 45, storm_at=storm)])
        runner = fake_runner_factory()
        for d in [f"2026-11-{x:02d}" for x in (12, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)]:
            ev, _, _ = run_day(base, d, gi, runner=runner)
        E = ev["2026-11-14_obs"]
        rep = evdir_file(base, E, "report_post7")
        check(E["status"] == "closed" and E["pre"]["status"] == "disabled"
              and not any(Path(a[1]).name == "dem_from_contours.py" for a in runner.calls)
              and rep is not None and "DEM_ENABLE=0" in rep.read_text()
              and "DEM_ENABLE=0" in evdir_file(base, E, "missed_storm").read_text(),
              "the cron's DEM_ENABLE=0: no DEM is built; the event runs its course and the emails say why")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "retry")
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(start, 45, storm_at=storm)])
        run_day(base, "2026-11-12", gi)
        tries = []
        st = []
        for d, rc in (("2026-11-15", 1), ("2026-11-16", 1), ("2026-11-17", 0)):
            ev, _, line = run_day(base, d, gi, mailer=lambda *a, rc=rc: tries.append(rc) or rc, email=True)
            m = ev["2026-11-14_obs"]["emails"][0]
            st.append((m["status"], m["tries"]))
        check(st == [("pending", 1), ("pending", 2), ("sent", 3)] and tries == [1, 1, 0],
              f"an email msmtp could not send stays pending and is retried on the next runs, then sent ({st})")
        texts = (line, (base / "logs/storm_watch.log").read_text(),
                 (base / "archive/storm_watch/state.json").read_text())
        check("missed_storm sent (2 recipients" in line and not any("@example.org" in t for t in texts),
              "the recipients are only counted: no address in the run line, logs/storm_watch.log or state.json")
        ev, _, _ = run_day(base, "2026-11-18", gi, mailer=lambda *a: tries.append(9) or 0, email=True)
        check(tries == [1, 1, 0], "a sent email is never sent again")
    # the first post-storm date (a September storm, frames from 09:30 UTC)
    E = {"observed": {"end": "2026-09-28T10:00:00Z", "onset": "2026-09-25T06:00:00Z"}, "forecast": None}
    ctx = type("C", (), {})()
    gi_s = {"latest": "2026-09-30", "earliest": {"2026-09-28": ep(datetime(2026, 9, 28, 9, 30, tzinfo=UTC))}}
    sw.set_first_post(ctx, E, gi_s)
    check(E["first_post"] == "2026-09-29" and E["first_post_exact"] and E["storm_dates"] == ["2026-09-25", "2026-09-28"],
          "a September end at 10Z with frames from 09:30Z: the end's date holds storm frames, first post date the next")
    gi_s["earliest"]["2026-09-28"] = ep(datetime(2026, 9, 28, 10, 30, tzinfo=UTC))
    sw.set_first_post(ctx, E, gi_s)
    check(E["first_post"] == "2026-09-28", "the same end, the first frame at 10:30Z: the end's date is a post date")
    # the wave record must cover a window before it is built
    cfg = sw.load_settings(None)[0]
    o = sw.Obs(sw.hour_of(datetime(2026, 11, 1, tzinfo=UTC)), 10 * 24)
    o.prim[:] = 0.8
    o.prim[3 * 24 + 2:3 * 24 + 10] = np.nan
    check(sw.wave_gaps(o, "2026-11-02", "2026-11-05", cfg) == ["2026-11-04"],
          "a date with fewer than 20 observed wave hours is a gap")
    Eg = {"waiting": {}}
    cg = type("C", (), {"cfg": cfg, "today": "2026-11-10"})()
    ok1 = sw.waves_ready(cg, Eg, "pre", o, "2026-11-02", "2026-11-05")
    cg.today = "2026-11-13"
    ok2 = sw.waves_ready(cg, Eg, "pre", o, "2026-11-02", "2026-11-05")
    check(not ok1[0] and ok2[0] and "gaps" in ok2[1], "a window with a wave gap waits, and is built with a flag after 3 days")


# --------------------------------------------------------------------------------------------
# T9b two storms close together, a forecast storm inside a compound event, camera-only waves, the state file,
#     an email that cannot be composed, the give-up rules
# --------------------------------------------------------------------------------------------

SHAPE = (2.0, 2.6, 3.0, 3.2, 3.1, 2.8, 2.6)


def wave_rows(start, n_days, hs_at):
    """[(epoch, hs, 'blend', qc 1, dir 60, tp 9)] hourly from start, hs_at(t) -> Hs."""
    return [(ep(t), hs_at(t), "blend", 1, 60, 9) for t in (start + timedelta(hours=k) for k in range(n_days * 24))]


def storms_at(*onsets, calm=0.9):
    """Hs(t): calm, with the SHAPE storm (then 20 h at 2.0 m) from each onset."""
    def hs(t):
        for s in onsets:
            k = (t - s).total_seconds() / 3600
            if 0 <= k < len(SHAPE):
                return SHAPE[int(k)]
            if len(SHAPE) <= k < len(SHAPE) + 20:
                return 2.0
        return calm
    return hs


def days_of(a, b):
    return sw.drange(a, b)


def body_of(base, E, kind):
    f = evdir_file(base, E, kind)
    return f.read_text() if f else ""


def t9b():
    print("T9b two storms close together, a compound event, camera-only waves, the state file, compose failures, "
          "give-up rules")
    start = datetime(2026, 10, 20, tzinfo=UTC)
    gi = simple_ground(days_of("2026-10-20", "2026-12-31"))
    # 1. a storm the forecast missed comes while the storm before has its post-window cut pending, and GNSS-R has
    #    stalled (the cut cannot be decided): it opens its own event at once
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "cutpend")
        A, B = datetime(2026, 11, 14, 3, tzinfo=UTC), datetime(2026, 11, 21, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv", wave_rows(start, 60, storms_at(A, B)))
        stalled = simple_ground(days_of("2026-10-20", "2026-11-19"))
        for d in days_of("2026-11-12", "2026-11-22"):
            ev, _, _ = run_day(base, d, stalled if d >= "2026-11-21" else gi, lag=0 if d >= "2026-11-21" else 2)
        Ea, Eb = ev["2026-11-14_obs"], ev.get("2026-11-21_obs")
        check(Ea["status"] == "post" and Ea["observed"]["cut_state"] == "pending" and Eb is not None
              and [m["kind"] for m in Eb["emails"]] == ["missed_storm"]
              and first_run_of(Eb, "OPENED") == "2026-11-21",
              f"a 3.2 m storm the forecast missed, inside the post window of the storm before (its cut pending: "
              f"GNSS-R stalled): its own event and email on the first run ({sorted(ev)})")
    # 2. a forecast storm while the storm before is still going on: its own watch event and heads-up at once
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "running")
        A0, A1 = datetime(2026, 11, 14, 3, tzinfo=UTC), datetime(2026, 11, 18, 0, tzinfo=UTC)
        B0, B1 = datetime(2026, 11, 20, 6, tzinfo=UTC), datetime(2026, 11, 21, 6, tzinfo=UTC)

        def hs_ab(t):
            if A0 <= t < A1:
                return 3.2 if t < A0 + timedelta(hours=12) else 2.0
            return 3.4 if B0 <= t < B1 else 0.9
        write_marconi(base / "archive/waves_marconi.csv", wave_rows(start, 40, hs_ab))
        pkB = datetime(2026, 11, 20, 18, tzinfo=UTC)
        for d in days_of("2026-11-13", "2026-11-17"):
            ev, _, line = run_day(base, d, gi, (lambda c: storm_fc(c, pkB, peak=5.0)) if d >= "2026-11-16" else None)
            if d == "2026-11-16":
                first = {k: [m["kind"] for m in e["emails"]] for k, e in ev.items()}
        fcs = [e for e in ev.values() if e["kind"] == "fc"]
        hb = body_of(base, fcs[0], "heads_up") if fcs else ""
        check(ev["2026-11-14_obs"]["status"] == "storm" and len(fcs) == 1 and first.get(fcs[0]["id"]) == ["heads_up"]
              and fcs[0]["prev_event"] == "2026-11-14_obs" and "whose storm is still going on" in hb,
              f"a storm forecast for 20 Nov while the 14 Nov storm is still rough: its own WATCH event with a heads-up "
              f"on the first triggering run (16 Nov), naming the storm before ({first})")
        check(len(fcs) == 1 and len(fcs[0]["forecast"]["runs"]) == 2
              and [m["kind"] for m in ev["2026-11-14_obs"]["emails"]] == ["missed_storm"],
              "... the next trigger updates that event, not the running storm")
        pp = (fcs[0]["pre_provisional"] or {}) if fcs else {}
        check(pp.get("status") == "none" and "no provisional pre-storm DEM: no date after the wave event" in hb
              and "provisional pre-storm DEM attached" not in hb,
              f"... and its heads-up has no provisional pre DEM of the beach before the running storm ({pp.get('why')})")
    # 3. a forecast storm whose sea comes as a wave event that the storm before merges before it confirms: closed as
    #    MERGED into the compound event (not a false alarm), with no pre DEM of its own
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "mergewatch")
        A, B0 = datetime(2026, 11, 14, 3, tzinfo=UTC), datetime(2026, 11, 17, 0, tzinfo=UTC)

        def hs_m(t):
            kb = (t - B0).total_seconds() / 3600
            if 0 <= kb < 72:
                return 2.1
            if 72 <= kb < 90:
                return 3.3
            if 90 <= kb < 100:
                return 2.0
            return storms_at(A)(t)
        write_marconi(base / "archive/waves_marconi.csv", wave_rows(start, 50, hs_m))
        runner = fake_runner_factory()
        pkB = datetime(2026, 11, 20, 9, tzinfo=UTC)
        for d in days_of("2026-11-12", "2026-11-26"):
            ev, _, _ = run_day(base, d, gi, (lambda c: storm_fc(c, pkB, peak=5.0)) if "2026-11-16" <= d <= "2026-11-19"
                               else None, runner=runner)
        Eb = ev["2026-11-18_fc"]
        Ea = ev["2026-11-14_obs"]
        mb = body_of(base, Eb, "merged")
        check(Eb["closed_reason"] == "merged" and [m["kind"] for m in Eb["emails"]] == ["heads_up", "merged"]
              and "part of the storm that event 2026-11-14_obs follows" in mb and "3.30 m" in mb
              and not Eb["pending"] and Eb["pre"]["status"] == "planned"
              and not any(f"{Eb['id']}/pre/" in a[3] for a in runner.calls if Path(a[1]).name == "dem_from_contours.py"),
              f"a forecast storm whose waves the storm before merged (a compound event): closed MERGED, not a false "
              f"alarm, with a short email and no pre DEM of its own ({Eb['closed_reason']})")
        check(any("COMPOUND" in n and Eb["id"] in n for n in Ea["notes"]) and Ea["storm_dates"] == ["2026-11-14", "2026-11-20"],
              "... the compound event's notes (its report emails) name the forecast event")
    # 4. a forecast storm a week after another: the provisional pre DEM never reaches back over the storm before
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "twoweek")
        A, B = datetime(2026, 11, 14, 3, tzinfo=UTC), datetime(2026, 11, 21, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv", wave_rows(start, 60, storms_at(A, B)))
        for d in days_of("2026-11-10", "2026-11-17"):
            dt = datetime.fromisoformat(d + "T12:50:00+00:00")
            fc = None
            for S in (A, B):
                if timedelta(0) < S - dt <= timedelta(days=4):
                    fc = (lambda S: (lambda c: storm_fc(c, S + timedelta(hours=6), peak=4.5)))(S)
            ev, _, _ = run_day(base, d, gi, fc)
        Eb = [e for e in ev.values() if e["opened_utc"][:10] == "2026-11-17"][0]
        pp = Eb["pre_provisional"] or {}
        hb = body_of(base, Eb, "heads_up")
        check(pp.get("status") == "none" and "no date after" in pp.get("why", "")
              and "no provisional pre-storm DEM" in hb and "2026-11-08 to 2026-11-14" not in hb,
              f"a heads-up 3 days after a storm (14 Nov), the ground file 2 days behind: no provisional pre DEM over "
              f"the storm before ({pp})")
    # 5. a short blip above the calm level days before a forecast storm does not shrink the planned pre window
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "blip")
        blip = datetime(2026, 11, 7, 15, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv",
                      wave_rows(start, 22, lambda t: 1.7 if blip <= t < blip + timedelta(hours=2) else 1.0))
        ev, _, _ = run_day(base, "2026-11-10", gi, lambda c: storm_fc(c, datetime(2026, 11, 13, 18, tzinfo=UTC)))
        E = list(ev.values())[0]
        check(E["pre"]["planned"] == ["2026-11-05", "2026-11-11"] and "pre-storm DEM: 2026-11-05 to 2026-11-11"
              in body_of(base, E, "heads_up"),
              f"a 2-hour 1.7 m blip on 7 Nov: the planned pre window stays 5-11 Nov ({E['pre']['planned']})")
        E2 = {"id": "x", "forecast": {"onset_fc": "2026-11-13T00:00:00Z"}, "observed": {"onset": None, "segments": []}}
        c = ctx_for(base, datetime(2026, 11, 14, 12, 50, tzinfo=UTC))
        seg = {"start": sw.hour_of(datetime(2026, 11, 12, 6, tzinfo=UTC)), "last": sw.hour_of(datetime(2026, 11, 12, 20,
               tzinfo=UTC)), "end": None, "peak": 2.2, "raw_confirm": False, "raw_max": None, "rough_hours": 15}
        check(sw.planned_pre_fc(c, {"events": []}, E2, [seg]) == ["2026-11-06", "2026-11-11"],
              "a wave event (2.2 m) starting the day before the forecast onset ends the planned window before it")
    # 5b. a thin planned window (a false alarm's pre DEM) is extended back, but never past the storm before or into
    #     a wave event that cuts windows
    for case in ("storm before", "wave event"):
        with tempfile.TemporaryDirectory() as tmp:
            base = make_base(tmp, "extend")
            c = ctx_for(base, datetime(2026, 11, 25, 12, 50, tzinfo=UTC), runner=fake_runner_factory())
            state = {"events": []}
            h = lambda d, hr: sw.hour_of(datetime.fromisoformat(f"{d}T{hr:02d}:00:00+00:00"))   # noqa: E731
            seg = {"start": h("2026-11-13", 3), "last": h("2026-11-14", 20), "end": h("2026-11-15", 0),
                   "peak": 3.2 if case == "storm before" else 2.2, "raw_confirm": False, "raw_max": None,
                   "rough_hours": 40, "n_confirm": 5}
            if case == "storm before":
                A = sw.new_event(c, state, "obs", "2026-11-13")
                A.update(status="closed", closed_reason="complete", storm_dates=["2026-11-13", "2026-11-14"])
                A["observed"]["segments"] = [sw.span_rec(seg)]
            E = sw.new_event(c, state, "fc", "2026-11-20")
            E["forecast"] = {"onset_fc": "2026-11-20T00:00:00Z", "end_fc": "2026-11-21T00:00:00Z"}
            g = simple_ground(days_of("2026-11-01", "2026-11-24"))
            for d in ("2026-11-15", "2026-11-16", "2026-11-17"):
                g["counts"][d] = {"c1": 1, "c2": 1}          # not usable
            sw.try_pre(c, state, E, [seg], g, None, window=["2026-11-16", "2026-11-19"])
            check(E["pre"]["window"][0] >= "2026-11-15" and E["pre"]["status"] == "inadequate",
                  f"a thin false-alarm pre window (2 usable dates) is not extended back past the {case} on 13-14 Nov "
                  f"({E['pre']['window']}, {E['pre']['status']})")
    # 5c. the provisional pre DEM: started after the storm before; none with fewer than 3 usable dates after it
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "prov")
        runner = fake_runner_factory()
        c = ctx_for(base, datetime(2026, 11, 19, 12, 50, tzinfo=UTC), runner=runner)
        state = {"events": []}
        h = lambda d, hr: sw.hour_of(datetime.fromisoformat(f"{d}T{hr:02d}:00:00+00:00"))   # noqa: E731
        seg = {"start": h("2026-11-13", 3), "last": h("2026-11-14", 20), "end": h("2026-11-15", 0), "peak": 3.2,
               "raw_confirm": False, "raw_max": None, "rough_hours": 40, "n_confirm": 5}
        A = sw.new_event(c, state, "obs", "2026-11-13")
        A.update(status="post", storm_dates=["2026-11-13", "2026-11-14"])
        A["observed"]["segments"] = [sw.span_rec(seg)]
        E = sw.new_event(c, state, "fc", "2026-11-20")
        E["forecast"] = {"onset_fc": "2026-11-20T00:00:00Z", "end_fc": "2026-11-21T00:00:00Z"}
        E["pre"]["planned"] = ["2026-11-15", "2026-11-19"]
        g = simple_ground(days_of("2026-11-01", "2026-11-18"))       # complete to 17 Nov
        sw.try_pre(c, state, E, [seg], g, None, provisional=True)
        ok3 = (E["pre_provisional"] or {}).get("window")
        g["counts"]["2026-11-16"] = {"c1": 1, "c2": 1}
        E["pre_provisional"] = None
        sw.try_pre(c, state, E, [seg], g, None, provisional=True)
        pp = E["pre_provisional"] or {}
        check(ok3 == ["2026-11-15", "2026-11-17"] and pp.get("status") == "none"
              and "2 usable date(s) in 2026-11-15..2026-11-17" in pp.get("why", ""),
              f"provisional pre DEM: the 7 dates to the last complete date, started after the storm before "
              f"({ok3}); none with only 2 usable dates after it ({pp.get('why')})")
    # 6. both buoys out, the camera alone sees the storm: closed UNVERIFIED, never a false alarm
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "camonly")
        s0, s1 = datetime(2026, 9, 13, 0, tzinfo=UTC), datetime(2026, 9, 15, 0, tzinfo=UTC)
        out0 = datetime(2026, 9, 10, tzinfo=UTC)
        rows = []
        for k in range(45 * 24):
            t = datetime(2026, 8, 20, tzinfo=UTC) + timedelta(hours=k)
            hs = 3.2 if s0 <= t < s1 else 0.9
            if t < out0:
                rows.append((ep(t), hs, "blend", 1, 60, 9))
            elif 10 <= t.hour <= 23:
                rows.append((ep(t), hs, "camera", 1, 60, 9))
        write_marconi(base / "archive/waves_marconi.csv", rows)
        gis = simple_ground(days_of("2026-08-20", "2026-10-31"))
        run_day(base, "2026-09-11", gis, lambda c: storm_fc(c, datetime(2026, 9, 14, tzinfo=UTC)))
        for d in days_of("2026-09-12", "2026-09-23"):
            ev, _, _ = run_day(base, d, gis)
        E = list(ev.values())[0]
        check(E["closed_reason"] == "unverified" and "false_alarm" not in [m["kind"] for m in E["emails"]],
              f"no buoy, camera-only rows at 3.2 m during the forecast storm: closed UNVERIFIED, not a false alarm "
              f"({E['closed_reason']})")
    # 7. a ground file with no georectified rows (or a replay without one): a watch event runs, no error
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "noground")
        write_marconi(base / "archive/waves_marconi.csv", wave_rows(start, 22, lambda t: 1.0))
        with contextlib.redirect_stdout(io.StringIO()):
            ev, _, line = run_day(base, "2026-11-10", {"latest": None, "counts": {}, "earliest": {}},
                                  lambda c: storm_fc(c, datetime(2026, 11, 13, 18, tzinfo=UTC)))
        E = list(ev.values())[0]
        check("ERROR" not in line and (E["pre_provisional"] or {}).get("status") == "none"
              and "no georectified rows" in E["pre_provisional"]["why"],
              "a ground file with no georectified rows: no provisional pre DEM, said why, no ERROR")
    # 8. --status is read-only; a missing or unreadable state file is restored from its previous copy
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "status")
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv", wave_rows(start, 45, storms_at(storm)))
        run_day(base, "2026-11-12", gi)
        run_day(base, "2026-11-15", gi)
        sdir = base / "archive/storm_watch"
        (sdir / "state.json").write_text('{"version": 1, "events": [')
        before = tree(base)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = sw.main(["--base", str(base), "--status"])
        check(rc == 0 and tree(base) == before and "unreadable" in out.getvalue() and "2026-11-14_obs" in out.getvalue(),
              "--status with an unreadable state file: says so, shows the events from the previous copy, moves nothing")
        ev, _, line = run_day(base, "2026-11-16", gi)
        check(list(ev) == ["2026-11-14_obs"] and any(p.name.startswith("state.json.bad-") for p in sdir.iterdir()),
              "the next daily run moves it aside and goes on with the restored events")
        (sdir / "state.json").unlink()
        ev, _, line = run_day(base, "2026-11-17", gi)
        check(list(ev) == ["2026-11-14_obs"] and "state file missing" in line,
              "state.json missing while state.json.prev exists: restored (never a first install that forgets the events)")
    # 9. an email that cannot be composed: a stand-in once, the real one on the next run
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "compose")
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv", wave_rows(start, 45, storms_at(storm)))
        run_day(base, "2026-11-12", gi)
        sent = []
        real = sw.compose

        def broken(ctx, state, E, m):
            if m["kind"] == "missed_storm":
                raise KeyError("peak_time")
            return real(ctx, state, E, m)
        sw.compose = broken
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                ev, _, line = run_day(base, "2026-11-15", gi, mailer=lambda s, t, b, a: sent.append((s, b)) or 0,
                                      email=True)
                m1 = dict(ev["2026-11-14_obs"]["emails"][0])
                ev, _, line2 = run_day(base, "2026-11-16", gi, mailer=lambda s, t, b, a: sent.append((s, b)) or 0,
                                       email=True)
        finally:
            sw.compose = real
        stamp = (base / "logs/storm_watch_last_run").read_text()
        ev, _, line3 = run_day(base, "2026-11-17", gi, mailer=lambda s, t, b, a: sent.append((s, b)) or 0, email=True)
        m = ev["2026-11-14_obs"]["emails"][0]
        check("ERROR: the missed_storm email of 2026-11-14_obs could not be composed (KeyError; try 1 of 3)" in line
              and m1["status"] == "pending" and m1.get("stub_sent") == "sent" and "try 2 of 3" in line2
              and stamp.startswith("error in: 2026-11-14_obs") and sent and "could not be composed" in sent[0][0],
              "an email that cannot be composed stays pending (not 'sent'); the run line and stamp say ERROR")
        check(len(sent) == 2 and len([s for s, _ in sent if "could not be composed" in s]) == 1 and m["status"] == "sent"
              and "A storm the forecast missed" in sent[-1][0] and "ERROR" not in line3,
              f"an email that cannot be composed: ERROR in the run line and stamp, a stand-in sent once, the real email "
              f"sent on the first run that composes it ({[s[:60] for s, _ in sent]})")
    # 10. the give-up rules
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "neverends")
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv",
                      wave_rows(start, 50, lambda t: 2.0 if t >= storm + timedelta(hours=7) else storms_at(storm)(t)))
        got = {}
        for d in ("2026-11-12", "2026-11-15", "2026-11-20", "2026-11-21", "2026-11-27", "2026-11-28"):
            ev, _, _ = run_day(base, d, gi)
            E = ev["2026-11-14_obs"] if "2026-11-14_obs" in ev else None
            got[d] = (E["status"], [m["kind"] for m in E["emails"]]) if E else None
        check(got["2026-11-20"] == ("storm", ["missed_storm"]) and got["2026-11-21"] == ("storm", ["missed_storm", "not_over"])
              and got["2026-11-27"][0] == "storm" and got["2026-11-28"] == ("closed", ["missed_storm", "not_over", "gave_up"])
              and E["closed_reason"] == "gave_up",
              f"a storm that never calms (peak 14 Nov 06Z): a not-over email 7 days after the peak, given up 14 days "
              f"after it ({got})")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "frozen")
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv", wave_rows(start, 50, storms_at(storm)))
        frozen = simple_ground(days_of("2026-10-20", "2026-11-13"))
        got = {}
        for d in ("2026-11-12", "2026-11-15", "2026-11-20", "2026-11-21", "2026-12-04", "2026-12-05"):
            ev, _, _ = run_day(base, d, frozen, lag=0)
            E = ev.get("2026-11-14_obs")
            got[d] = (E["status"], [m["kind"] for m in E["emails"]][1:]) if E else None
        wb, gb = body_of(base, E, "waiting_ground"), body_of(base, E, "gave_up")
        check(got["2026-11-20"] == ("post", []) and got["2026-11-21"] == ("post", ["waiting_ground"])
              and got["2026-12-04"] == ("post", ["waiting_ground"]) and got["2026-12-05"] == ("closed", ["waiting_ground", "gave_up"])
              and E["closed_reason"] == "gave_up",
              f"the ground file frozen at 13 Nov (the pre needs 13 Nov complete): a waiting email 8 days on, given up "
              f"after STORM_MAX_WAIT_DAYS (21) ({got})")
        check("needs the waterline points of 2026-11-13, but the ground file's last complete date is 2026-11-12" in wb
              and "last complete date is 2026-11-12" in gb and "is closed: gave up waiting" in gb,
              "... the emails name the ground file's last complete date (12 Nov), not the date it waits for")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "receding")
        write_marconi(base / "archive/waves_marconi.csv", wave_rows(start, 45, lambda t: 1.0))
        got = []
        for d in days_of("2026-11-10", "2026-11-22"):
            ev, _, _ = run_day(base, d, gi, lambda c: storm_fc(c, c.now + timedelta(days=3, hours=12)))
            E = sorted(ev.values(), key=lambda e: e["opened_utc"])[0]
            got.append((d, E["status"]))
        closed = [d for d, s in got if s == "closed"]
        check(closed and closed[0] == "2026-11-21" and E["closed_reason"] == "false_alarm"
              and "kept receding" in body_of(base, E, "false_alarm"),
              f"a forecast onset that recedes every day: closed as a false alarm on day 11 ({closed[:1]})")


def first_run_of(E, prefix):
    return next((h[:10] for h in E["history"] if h[12:].startswith(prefix)), None)


# --------------------------------------------------------------------------------------------
# T10 email
# --------------------------------------------------------------------------------------------

def fake_msmtp(tmp, rc=0):
    b = Path(tmp) / "mbin"
    b.mkdir(exist_ok=True)
    (b / "msmtp").write_text(f'#!/bin/bash\necho "$@" >> "{tmp}/msmtp.argv"\ncat > "{tmp}/msmtp.stdin"\n'
                             f'[ {rc} -ne 0 ] && echo "msmtp: recipient address $1 not accepted by the server" >&2\n'
                             f'exit {rc}\n')
    (b / "msmtp").chmod(0o755)
    return b


def fixture_routes():
    """FakeNet routes serving the SYNTHETIC fixtures for every source a live run fetches."""
    text = json.dumps({"properties": {"periods": [{"name": "Tonight", "detailedForecast": "Northeast wind 25 to 30 "
                                                                                          "kt. Seas 8 to 11 ft."}]}})
    return [("/points/", [(200, (FIX / "nws_points_SYNTHETIC.json").read_bytes())]),
            ("api.weather.gov/gridpoints/BOX/101,101/forecast", [(200, text.encode())]),
            ("/gridpoints/", [(200, (FIX / "nws_grid_SYNTHETIC.json").read_bytes())]),
            ("/alerts/", [(200, (FIX / "nws_alerts_SYNTHETIC.json").read_bytes())]),
            ("/zones/", [(200, b"{}")]),
            ("marine-api", [(200, (FIX / "om_marine_SYNTHETIC.json").read_bytes())]),
            ("open-meteo", [(200, (FIX / "om_wind_SYNTHETIC.json").read_bytes())]),
            ("polar.ncep", [(200, (FIX / "gfswave_44008_SYNTHETIC.bull").read_bytes())])]


def t10():
    print("T10 email")
    with tempfile.TemporaryDirectory() as tmp:
        b = fake_msmtp(tmp)
        png = Path(tmp) / "map.png"
        png.write_bytes(b"\x89PNG\r\n\x1a\nfake")
        tif = Path(tmp) / "dem.tif"
        tif.write_bytes(b"II*\x00fake")
        old = os.environ["PATH"]
        os.environ["PATH"] = f"{b}:{old}"
        try:
            rc = sw.send("[CACO-05 storm watch] test", ["a@example.org", "b@example.org"], "body\n", [png, tif], False,
                         waits=(0, 0, 0))
        finally:
            os.environ["PATH"] = old
        argv = (Path(tmp) / "msmtp.argv").read_text().split()
        msg = (Path(tmp) / "msmtp.stdin").read_text(errors="replace")
        check(rc == 0 and argv == ["a@example.org", "b@example.org"], "msmtp gets each recipient as its own argument")
        check("Subject: [CACO-05 storm watch] test" in msg and "Content-Type: image/png" in msg
              and "Content-Type: image/tiff" in msg and "From:" not in msg,
              "subject prefix; a PNG as image/png, a GeoTIFF as image/tiff; no From (msmtp's account sets it)")
        (Path(tmp) / "x").mkdir()
        b2 = fake_msmtp(Path(tmp) / "x", rc=1)
        os.environ["PATH"] = f"{b2}:{old}"
        err = io.StringIO()
        try:
            with contextlib.redirect_stderr(err):
                rc = sw.send("s", ["a@example.org"], "b", [], False, waits=(0, 0, 0))
        finally:
            os.environ["PATH"] = old
        check(rc == 1, "a failing msmtp, tried 3 times: rc 1")
        check("not accepted" in err.getvalue() and "<address>" in err.getvalue() and "@example.org" not in err.getvalue(),
              "msmtp's error text is printed (into the detail log) with every address blanked out")
    # --dry-run writes nothing and sends nothing
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(datetime(2026, 11, 1, tzinfo=UTC), 20)])
        fdir = Path(tmp) / "fc"
        fdir.mkdir()
        for name, src in (("nws_grid.json", "nws_grid_SYNTHETIC.json"), ("nws_alerts.json", "nws_alerts_SYNTHETIC.json"),
                          ("om_marine.json", "om_marine_SYNTHETIC.json"), ("om_wind.json", "om_wind_SYNTHETIC.json"),
                          ("gfswave_44008.bull", "gfswave_44008_SYNTHETIC.bull"),
                          ("nws_points.json", "nws_points_SYNTHETIC.json")):
            shutil.copy(str(FIX / src), str(fdir / name))
        b = fake_msmtp(tmp)
        before = tree(base)
        old = os.environ["PATH"]
        os.environ["PATH"] = f"{b}:{old}"
        out = io.StringIO()
        try:
            with contextlib.redirect_stdout(out):
                rc = sw.main(["--base", str(base), "--dry-run", "--forecast-dir", str(fdir), "--now",
                              "2026-11-16T12:50:00Z"])
        finally:
            os.environ["PATH"] = old
        text = out.getvalue()
        check(rc == 0 and tree(base) == before and not (Path(tmp) / "msmtp.argv").exists(),
              "--dry-run: the base folder is byte-for-byte as before, and msmtp is never called")
        check("TRIGGER" in text and "would send to" in text and "heads_up" in text and "DRY RUN" in text,
              "--dry-run prints the decision and the heads-up it would send")
        with contextlib.redirect_stdout(io.StringIO()):
            rc1 = sw.main(["--base", str(base), "--now", "2026-11-16T12:50:00Z"])
            rc2 = sw.main(["--base", str(base), "--dry-run", "--now", "2026-11-16 12:50"])
        check(rc1 == 1 and rc2 == 1, "--now only with --dry-run, and only with Z or an offset")
    # a dry run that FETCHES (the guide's first command on the NUC, ./storm_watch.sh --dry-run --save-raw DIR)
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(datetime(2026, 11, 1, tzinfo=UTC), 20)])
        raw = Path(tmp) / "sw_raw"
        net = FakeNet(fixture_routes(), default=(404, b"{}"))
        before = tree(base)
        ctx = ctx_for(base, NOW, mode="dry", fetch=True, save_raw=str(raw), opener=net)
        with contextlib.redirect_stdout(io.StringIO()):
            rc, line = sw.run_once(ctx)
        check(rc == 0 and len(net.seen) >= 7 and "TRIGGER" in line and tree(base) == before,
              "--dry-run that fetches every source: the base folder is byte-for-byte as before (no forecast archive)")
        saved = sorted(p.name for p in raw.iterdir()) if raw.exists() else []
        check(set(sw.SOURCE_FILES) <= set(saved), f"--save-raw keeps each response, plain ({saved})")
        raw.mkdir(exist_ok=True)
        for p in raw.iterdir():
            os.utime(p, (NOW.timestamp(), NOW.timestamp()))      # as if saved at NOW (real_captures reads the time)
        n = len(FAILED)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            real_captures(str(raw))
        failed_here = FAILED[n:]
        del FAILED[n:]
        check(saved and not failed_here and "trigger: True" in out.getvalue(),
              f"... and those captures pass test_storm_watch.py --real-captures ({failed_here})")


# --------------------------------------------------------------------------------------------
# T11 replay
# --------------------------------------------------------------------------------------------

def t11():
    print("T11 replay")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        (base / "archive/gnssr_record_start.txt").write_text("2026-08-01 90\n")
        storm = datetime(2026, 9, 23, 3, tzinfo=UTC)
        rows = storm_hours(datetime(2026, 8, 20, tzinfo=UTC), 70, storm_at=storm,
                           storm_hs=(2.4, 3.2, 3.8, 4.2, 4.6, 4.4, 4.0, 3.8, 3.6), calm=1.0, tail_hours=40, tail_hs=2.6)
        write_buoy(base / "archive/waves_44008.csv", [(e, h, 10.0, 60) for e, h in rows])
        days = [sw.dadd("2026-09-01", k) for k in range(50)]
        by_cam = ground_rows_simple(days, 6, hs=lambda d, k: 0.8)
        write_ground_file(base / "contour_points_ground.csv", by_cam)
        b = fake_msmtp(tmp)
        out_dir = Path(tmp) / "replay_out"
        cron, cfg = cfg_for(base)
        runner = fake_runner_factory()
        before = tree(base)
        old = os.environ["PATH"]
        os.environ["PATH"] = f"{b}:{old}"
        try:
            with contextlib.redirect_stdout(io.StringIO()) as so:
                rc = sw.replay(base, "2026-09-15", "2026-10-12", out_dir, cfg, cron, runner=runner, local_now=LOCAL,
                               sleep=nosleep)
        finally:
            os.environ["PATH"] = old
        printed = so.getvalue()
        check(rc == 0 and tree(base) == before, "a replay writes nothing under the station folder")
        st = json.loads((out_dir / "state.json").read_text())
        ev = st["events"]
        check(len(ev) == 1 and ev[0]["kind"] == "fc" and ev[0]["hindcast"]
              and ev[0]["forecast"]["sources"].get("hindcast", "").startswith("observed waves"),
              f"the storm of 23 Sep opened from the perfect forecast (observed waves), as a hindcast ({[e['id'] for e in ev]})")
        if ev:
            opened = ev[0]["opened_utc"][:10]
            lead = (storm - datetime.fromisoformat(opened + "T13:50:00+00:00")).total_seconds() / 86400
            check(3.0 <= lead <= 4.2, f"opened {lead:.1f} days ahead of the onset (the 96 h lead)")
            check(ev[0]["status"] == "closed" and ev[0]["closed_reason"] == "complete"
                  and ev[0]["pre"]["window"][1] == "2026-09-22" and ev[0]["post7"],
                  f"the event runs its course: pre ending 22 Sep, post7 {ev[0]['post7'].get('window') if ev[0]['post7'] else None}")
            check(all(m["status"] == "dry" for m in ev[0]["emails"]) and not (Path(tmp) / "msmtp.argv").exists()
                  and all("[HINDCAST]" in m["subject"] for m in ev[0]["emails"]),
                  "replay emails: bodies only, subject [HINDCAST], msmtp never called")
        check((out_dir / "replay_summary.json").exists() and "REPLAY 2026-09-15 to 2026-10-12" in printed
              and (out_dir / "storm_watch.log").exists() and not (base / "logs/storm_watch.log").exists(),
              "replay_summary.json, its own log and a printed table")
        r = lambda *a, **k: sw.replay_refusal(base, *a, **k)   # noqa: E731
        check("inside" in r("2026-09-15", "2026-09-20", base / "archive/storm_watch/x", cron, cfg, 2),
              "refused: an --out-dir inside archive/storm_watch")
        check("inside" in r("2026-09-15", "2026-09-20", base / "archive/dems/r", cron, cfg, 2)
              and "station folder" in r("2026-09-15", "2026-09-20", base, cron, cfg, 2),
              "refused: an --out-dir inside archive/dems, or the station folder itself")
        check("state.json exists" in r("2026-09-15", "2026-09-20", out_dir, cron, cfg, 2),
              "refused: an --out-dir that already holds a replay")
        check("before the live GNSS-R record" in r("2026-07-15", "2026-07-20", Path(tmp) / "o2", cron, cfg, 2),
              "refused: a START before the GNSS-R record start")
        (base / "archive/gnssr_record_start.txt").write_text("2024-01-01 900\n")
        msg = r("2024-11-05", "2024-11-20", Path(tmp) / "o3", cron, cfg, 2)
        check("2024-10-23" in msg and "blind" in msg,
              "refused: a range whose padding reaches 2024-10-23 (the blind check's date), never read")
        (base / "archive/gnssr_record_start.txt").unlink()
        check("no GNSS-R record start" in r("2026-09-15", "2026-09-20", Path(tmp) / "o4", cron, cfg, 2),
              "refused: no record-start file")
        t = sw.live_run_time("2026-07-01"), sw.live_run_time("2026-12-01")
        check(t[0].hour == 12 and t[1].hour == 13 and t[0].minute == 50, "the live run time: 12:50Z in CDT, 13:50Z in CST")
    # the forecast archive: a live run keeps each response gzipped, prunes after 730 days; a replay reads it back
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(datetime(2026, 10, 1, tzinfo=UTC), 60, calm=1.0)])
        fdir = base / "archive/storm_watch/forecasts"
        old_d, keep_d = fdir / sw.dadd("2026-11-16", -731), fdir / sw.dadd("2026-11-16", -730)
        for d in (old_d, keep_d):
            d.mkdir(parents=True)
            (d / "1250Z_om_marine.json.gz").write_bytes(gzip.compress(b"{}"))
        ctx = ctx_for(base, NOW, fetch=True, opener=FakeNet(fixture_routes(), default=(404, b"{}")),
                      runner=fake_runner_factory())
        with contextlib.redirect_stdout(io.StringIO()):
            rc, line = sw.run_once(ctx)
        today = fdir / "2026-11-16"
        arch = sorted(p.name for p in today.iterdir()) if today.exists() else []
        bull = today / "1250Z_gfswave_44008.bull.gz"
        check(rc == 0 and "TRIGGER" in line and {f"1250Z_{n}.gz" for n in sw.SOURCE_FILES} <= set(arch)
              and bull.exists() and gzip.decompress(bull.read_bytes()) == (FIX / "gfswave_44008_SYNTHETIC.bull").read_bytes(),
              f"a live run archives each response gzipped in forecasts/<date>/<HHMMZ>_<name>.gz ({len(arch)} files)")
        check(not old_d.exists() and keep_d.exists() and today.exists(),
              "forecasts older than STORM_FORECAST_KEEP_DAYS (730) are pruned; 730 days old and today's are kept")
        (base / "archive/gnssr_record_start.txt").write_text("2026-08-01 90\n")
        d17 = fdir / "2026-11-17"
        d17.mkdir(parents=True)
        (d17 / "1250Z_gfswave_44008.bull.gz").write_bytes(gzip.compress((FIX / "gfswave_44008_SYNTHETIC.bull")
                                                                        .read_bytes()))   # only the bulletin that day
        cron, cfg = cfg_for(base)
        out = Path(tmp) / "rp"
        with contextlib.redirect_stdout(io.StringIO()):
            rc = sw.replay(base, "2026-11-16", "2026-11-18", out, cfg, cron, runner=fake_runner_factory(),
                           local_now=LOCAL, sleep=nosleep)
        log = (out / "storm_watch.log").read_text().splitlines()
        st = json.loads((out / "state.json").read_text())
        check(rc == 0 and len(log) == 3 and log[0].startswith("[2026-11-16T13:50:00Z]") and "om_marine ok" in log[0]
              and "hindcast" not in log[0] and "TRIGGER" in log[0] and st["events"]
              and st["events"][0]["opened_utc"].startswith("2026-11-16"),
              f"a replay day with archived forecasts uses them, not the observed waves: {log[0][23:120] if log else ''}")
        check(len(log) == 3 and "gfswave_44008" in log[1] and "hindcast" not in log[1] and "forecast hindcast" in log[2],
              "a day with only the GFS-Wave bulletin archived reads it; a day with nothing archived: the perfect forecast")
        check(not any("ERROR" in ln for ln in log),
              "a replay without a ground file follows its events with no ERROR (no DEM is built)")


# --------------------------------------------------------------------------------------------
# T12 climatology
# --------------------------------------------------------------------------------------------

def t12():
    print("T12 climatology")
    cfg = sw.load_settings(None)[0]
    lines = ["#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS  TIDE",
             "#yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC   mi    ft"]
    t = datetime(2025, 1, 1, tzinfo=UTC)
    for k in range(24 * 30):
        tt = t + timedelta(hours=k)
        hs, d = 1.0, 60
        if 100 <= k < 107:
            hs, d = 3.8, 60
        elif 300 <= k < 308:
            hs, d = 4.2, 230
        elif 500 <= k < 502:
            hs, d = 3.6, 60
        lines.append(f"{tt:%Y %m %d %H} 00 MM MM MM {hs:5.2f} 10.00 7.00 {d:3d} 1010.0 MM MM MM MM MM")
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "44008h2025.txt.gz"
        p.write_bytes(gzip.compress("\n".join(lines).encode()))
        with contextlib.redirect_stdout(io.StringIO()):
            rc, res = sw.climatology(cfg, [str(p)])
        n35, n40 = res["thresholds"]["3.5"]["events"], res["thresholds"]["4.0"]["events"]
        check(rc == 0 and n35 == 1 and n40 == 0,
              f"synthetic stdmet: 7 h of 3.8 m from 60 counts; 8 h of 4.2 m from 230 and 2 h of 3.6 m do not "
              f"(3.5 m: {n35}, 4.0 m: {n40})")
    with contextlib.redirect_stdout(io.StringIO()):
        rc, res = sw.climatology(cfg, [str(HERE / "sig1000_waves_ALL.csv")])
    c = {k: len(v["catches"]) for k, v in res["thresholds"].items()}
    check(len(res["opening_segments"]) == 5 and c["3.5"] == 3 and c["3.0"] == 4 and c["4.0"] == 0,
          f"the real ADCP / 0.83 as a perfect forecast, >= 3 h: 3.5 m catches 3 of the 5 storms that open an event, "
          f"3.0 m 4 (finding 7: {c})")
    cfg6 = dict(cfg, STORM_MIN_HOURS=6)
    cfg1 = dict(cfg, STORM_MIN_HOURS=1)
    with contextlib.redirect_stdout(io.StringIO()):
        r6 = sw.climatology(cfg6, [str(HERE / "sig1000_waves_ALL.csv")])[1]
        r1 = sw.climatology(cfg1, [str(HERE / "sig1000_waves_ALL.csv")])[1]
    check(len(r6["thresholds"]["3.5"]["catches"]) == 1 and len(r1["thresholds"]["3.5"]["catches"]) == 5,
          "3.5 m for >= 6 h catches 1 of 5, for >= 1 h 5 of 5 (why STORM_MIN_HOURS is 3)")


# --------------------------------------------------------------------------------------------
# T13 storm_watch.sh
# --------------------------------------------------------------------------------------------

def t13():
    print("T13 storm_watch.sh")
    r = subprocess.run(["bash", "-n", str(HERE / "storm_watch.sh")], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    check(r.returncode == 0, "bash -n storm_watch.sh")
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp) / "station"
        base.mkdir()
        shutil.copy(str(HERE / "storm_watch.sh"), str(base))
        b = Path(tmp) / "bin"
        b.mkdir()
        tab = Path(tmp) / "crontab.txt"
        tab.write_text("30 12 * * * /mnt/I2Rgus_Data/waterline/waterline_timex_cron.sh\n"
                       "40 6-20 * * * /mnt/x/owg.sh # owg.sh (optical wave gauge): hourly update\n")
        (b / "crontab").write_text(f'#!/bin/bash\nif [ "$1" = "-l" ]; then cat "{tab}"; else cat > "{tab}"; fi\n')
        (b / "crontab").chmod(0o755)
        env = dict(os.environ, PATH=f"{b}:{os.environ['PATH']}")
        sh = str(base / "storm_watch.sh")
        r1 = subprocess.run([sh, "--install-cron"], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
        r2 = subprocess.run([sh, "--install-cron"], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
        lines = tab.read_text().splitlines()
        sw_lines = [ln for ln in lines if "storm_watch.sh" in ln]
        check(r1.returncode == 0 and r2.returncode == 0 and len(sw_lines) == 1
              and sw_lines[0].startswith(f"50 7 * * * {base}/storm_watch.sh # storm_watch.sh (storm watch)")
              and len(lines) == 3, f"--install-cron adds the 07:50 line once (twice run: still one); others kept ({lines})")
        saved = (base / "logs/crontab.before_storm_watch_install").read_text()
        check("owg.sh" in saved and "waterline_timex_cron.sh" in saved, "the previous crontab saved in logs/")
        r3 = subprocess.run([sh, "--remove-cron"], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
        lines = tab.read_text().splitlines()
        check(r3.returncode == 0 and not any("storm_watch" in ln for ln in lines) and len(lines) == 2,
              "--remove-cron takes it out, and only it")
        tab.write_text("0 1 * * * /something/else.sh\n")
        r4 = subprocess.run([sh, "--install-cron"], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
        check(r4.returncode == 1 and tab.read_text() == "0 1 * * * /something/else.sh\n",
              "no waterline_timex_cron.sh line in the crontab (another user?): refused, unchanged")
        (b / "crontab").write_text("#!/bin/bash\nexit 1\n")
        r5 = subprocess.run([sh, "--install-cron"], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
        check(r5.returncode == 1 and "not changing" in r5.stdout, "crontab -l failing: refused (nothing wiped)")
        (b / "python3").write_text("#!/bin/bash\necho boom; exit 2\n")
        (b / "python3").chmod(0o755)
        r6 = subprocess.run([sh], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                            universal_newlines=True)
        runlog = (base / "logs/storm_watch.log").read_text()
        check(r6.returncode == 2 and "ERROR: storm_watch.py exited 2" in runlog
              and "boom" in (base / "logs/storm_watch_detail.log").read_text(),
              "a run that fails: the ERROR line in logs/storm_watch.log, its output in storm_watch_detail.log")
        (b / "python3").write_text("#!/bin/bash\necho '--- would send to: a@example.org b@example.org ---'; exit 0\n")
        shown = []
        for opt in ("--dry-run", "--status"):
            r7 = subprocess.run([sh, opt], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                stdin=subprocess.DEVNULL, universal_newlines=True)
            shown.append(r7.returncode == 0 and "a@example.org" in r7.stdout)
        detail = (base / "logs/storm_watch_detail.log").read_text()
        check(all(shown) and "@example.org" not in detail and "@example.org" not in (base / "logs/storm_watch.log").read_text(),
              "--dry-run / --status print the recipients on the screen only: never into storm_watch_detail.log")
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp)
        p = hold_lock(base / "logs/storm_watch.lock")
        try:
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                rc = sw.main(["--base", str(base), "--no-fetch", "--no-email"])
        finally:
            drop_lock(p)
        log = (base / "logs/storm_watch.log").read_text()
        check(rc == 0 and "SKIPPED: another storm watch run holds logs/storm_watch.lock" in log
              and not (base / "archive/storm_watch").exists(),
              "storm_watch.py's own lock held by another process: one SKIPPED line, exit 0, nothing else")


# --------------------------------------------------------------------------------------------
# T14 a station other than CACO-05: the defaults table, the settings files, --check-settings, the station's
#     DEM_MAX_HS, recipients, the buoys, the Open-Meteo facing test, the time zone
# --------------------------------------------------------------------------------------------

# A SYNTHETIC station's storm_watch_station.conf, as coastcam-station-setup would write it (a west-facing beach; NDBC
# ids and zones of that coast; no network is used)
STATION_TEXT = ('# storm_watch_station.conf -- TEST-01 (SYNTHETIC, for the tests)\n'
                'STORM_STATION=TEST-01\nSTORM_SITE="Test Beach"\nSTORM_TZ=America/Los_Angeles\n'
                'STORM_POINT=37.75,-122.51\nSTORM_MARINE_POINT=37.75,-122.70\nSTORM_BUOY_POINT=37.76,-122.83\n'
                'STORM_BUOY=46026\nSTORM_BUOY2=46237\nSTORM_ALERT_ZONES=PZZ545,CAZ006\nSTORM_FACING=270\n'
                'STORM_MIN_CELL_LON=""\nSTORM_DIR_SECTOR=190-350\n')


def station_base(tmp, name, extra="", settings=None):
    """A station folder with STATION_TEXT (+ extra) as its storm_watch_station.conf, owg.sh present."""
    base = make_base(tmp, name, settings=settings)
    (base / "storm_watch_station.conf").write_text(STATION_TEXT + extra)
    return base


def station_cfg(base):
    cron, _ = sw.cron_settings(base)
    return cron, sw.load_settings([base / "storm_watch_station.conf", base / "storm_watch_settings.conf"], cron)[0]


def check_settings_run(base):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = sw.main(["--base", str(base), "--check-settings"])
    return rc, out.getvalue()


def t14():
    print("T14 another station: defaults table, settings files, --check-settings, DEM_MAX_HS, recipients, buoys, "
          "facing, time zone")
    # 1. no buoy number outside the one defaults table (correction C10)
    src = (HERE / "storm_watch.py").read_text().splitlines()
    starts = [i for i, ln in enumerate(src) if ln.startswith("# >>> CACO-05 defaults")]
    ends = [i for i, ln in enumerate(src) if ln.startswith("# <<< CACO-05 defaults")]
    a, b = (starts[0], ends[0]) if len(starts) == 1 and len(ends) == 1 else (-1, -1)
    lit = re.compile(r"440(08|13)")
    outside = [i + 1 for i, ln in enumerate(src) if lit.search(ln) and not a < i < b]
    inside = "\n".join(src[a:b]) if a >= 0 else ""
    check(len(starts) == 1 and len(ends) == 1 and not outside and "44008" in inside and "44013" in inside,
          f"storm_watch.py names 44008 / 44013 only inside its one CACO-05 defaults table (lines {a + 1}-{b + 1}; "
          f"outside: {outside})")
    cfg0 = sw.load_settings(None)[0]
    check(sw.buoys(cfg0) == ("44008", "44013") and cfg0["STORM_GFSWAVE_STATION"] == "44008" and cfg0["STORM_MODE"] == "full"
          and cfg0["STORM_STATION"] == "CACO-05" and cfg0["STORM_FACING"] is None and cfg0["STORM_MIN_CELL_LON"] == -69.93
          and sw.USER_AGENT == "caco05-storm-watch (USGS coastal camera research)" and sw.SUBJECT == "[CACO-05 storm watch]"
          and sw.SOURCE_FILES[-1] == "gfswave_44008.bull" and sw.SRC_LABELS[2:4] == ["44013 converted", "44008 x 0.83"],
          "no station file: CACO-05's buoys, mode full, the User-Agent, subject, file names and labels of today")

    # 2. the settings files: storm_watch_station.conf, then storm_watch_settings.conf (the last value wins)
    with tempfile.TemporaryDirectory() as tmp:
        base = station_base(tmp, "prec", extra="STORM_MODE=alerts\nSTORM_HS_M=4.0\nSTORM_CONFIRM_HS=2.2\n",
                            settings='STORM_HS_M=4.5\nSTORM_SITE="Own Beach"\n')
        cron, cfg = station_cfg(base)
        check(cfg["STORM_HS_M"] == 4.5 and cfg["STORM_CONFIRM_HS"] == 2.2 and cfg["STORM_SITE"] == "Own Beach"
              and cfg["STORM_BUOY"] == "46026" and cfg["STORM_MODE"] == "alerts" and cfg["STORM_FACING"] == 270.0
              and cfg["STORM_MIN_CELL_LON"] is None and cfg["STORM_GFSWAVE_STATION"] == "46026",
              "storm_watch_station.conf then storm_watch_settings.conf: the settings file's value wins (STORM_HS_M "
              "4.5), the station file's others stand")
        rc, out = check_settings_run(base)
        check(rc == 0 and "site Own Beach" in out and "mode alerts" in out and "buoys 46026 and 46237" in out,
              "main reads both files in that order (--check-settings shows the site from the settings file, the mode "
              "and buoys from the station file)")
        (base / "storm_watch_settings.conf").write_text("STORM_MODE=full\n")
        with contextlib.redirect_stdout(io.StringIO()) as so:
            sw.main(["--base", str(base), "--dry-run", "--no-fetch", "--now", "2026-11-16T12:50:00Z"])
        (base / "storm_watch_settings.conf").write_text("")
        with contextlib.redirect_stdout(io.StringIO()) as so2:
            sw.main(["--base", str(base), "--dry-run", "--no-fetch", "--now", "2026-11-16T12:50:00Z"])
        check("alerts mode" not in so.getvalue() and "alerts mode (no DEM)" in so2.getvalue(),
              "a run takes STORM_MODE from the settings file over the station file's")

    # 3. --check-settings: exit 1 on an unknown key, a value that does not parse, 'KEY = v', a missing station fact
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "chk", settings=None)
        rc, out = check_settings_run(base)
        check(rc == 0 and "settings OK" in out and "WARN" not in out and "owg.sh RECIPIENTS" in out,
              "CACO-05 without any settings file: OK, no warning (recipients from owg.sh, only counted)")
        check("@" not in out, "--check-settings never prints an address")
        cases = [("STORM_FOO=1\n", "STORM_FOO is not a storm watch setting"),
                 ("STORM_HS_M=abc\n", "does not parse"),
                 ("STORM_HS_M = 4.0\n", "no spaces around '='"),
                 ("STORM_HS_M= 4.0\n", "no spaces around '='"),
                 ("STORM_MODE=partial\n", "STORM_MODE='partial' does not parse"),
                 ("STORM_TZ=Mars/Olympus_Mons\n", "STORM_TZ='Mars/Olympus_Mons' does not parse"),
                 ("STORM_BUOY=44 008\n", "text after the value"),
                 ('STORM_SITE="Unclosed\n', "quote is not closed")]
        got = []
        for text, want in cases:
            (base / "storm_watch_settings.conf").write_text("# a comment\n\n" + text)
            rc, out = check_settings_run(base)
            got.append((text.strip(), rc, want in out and "settings REFUSED" in out))
        check(all(rc == 1 and ok for _, rc, ok in got), f"--check-settings exit 1: unknown key, unparsable values, "
                                                          f"'KEY = v', 'KEY= v', text after a value, an open quote ({got})")
        (base / "storm_watch_settings.conf").write_text("STORM_STATION=TEST-02\n")
        rc, out = check_settings_run(base)
        warned = [k for k in sw.MARCONI_FACTS if f"WARN: {k} is CACO-05's (Marconi Beach) value" in out]
        check(rc == 0 and set(warned) == set(sw.MARCONI_FACTS) and "thresholds" in out and "--emit-settings" in out,
              f"another station on Marconi's facts: exit 0, a WARN for each ({len(warned)} of "
              f"{len(sw.MARCONI_FACTS)}) and for the thresholds")
        (base / "storm_watch_settings.conf").unlink()
        (base / "storm_watch_station.conf").write_text(STATION_TEXT.replace("STORM_FACING=270\n", ""))
        rc, out = check_settings_run(base)
        check(rc == 1 and "the station fact STORM_FACING is missing" in out,
              "a storm_watch_station.conf without one of its station facts (STORM_FACING): exit 1")
        (base / "storm_watch_station.conf").write_text(STATION_TEXT)
        rc, out = check_settings_run(base)
        check(rc == 0 and "CACO-05's (Marconi Beach) value" not in out,
              "a complete station file with its own facts: exit 0, no Marconi-fact warning")
        check("emails ON: 0 recipient(s)" in out and "no owg.sh fallback" in out and "WARN: emails are on" in out,
              "... STORM_EMAIL=1 without STORM_EMAIL_TO on a setup station: 0 recipients (no owg.sh), a WARN")
        (base / "storm_watch_station.conf").write_text(STATION_TEXT.replace("STORM_BUOY2=46237", "STORM_BUOY2=44013"))
        rc, out = check_settings_run(base)
        check(rc == 0 and "WARN: STORM_BUOY2 is CACO-05's (Marconi Beach) value '44013'" in out,
              "a station file giving one of Marconi's facts (STORM_BUOY2 44013): a WARN, exit 0")
        (base / "storm_watch_station.conf").write_text(STATION_TEXT)
        (base / "storm_watch_settings.conf").write_text("STORM_HS_M=4.5\nSTORM_EMAIL=0\n")
        cron, cfg = station_cfg(base)
        rc, out = check_settings_run(base)
        check(rc == 0 and cfg["STORM_HS_M"] == 4.5 and "emails off" in out, "station file and settings file together: OK")

    # 4. waterline_station.conf's DEM_MAX_HS is applied after the cron's lines
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "overlay")
        c0, _ = sw.cron_settings(base)
        (base / "waterline_station.conf").write_text('WATERLINE_STAGE=""\nOWG_ENABLE="0"\nDEM_MAX_HS="1.0"\n'
                                                     'DEM_CELL="9.0"\nSTATION_NAME="TEST-01"\nexport STATION_NAME\n')
        cron, notes = sw.cron_settings(base)
        cfg = sw.load_settings(None, cron)[0]
        check(c0["DEM_MAX_HS"] == "1.5" and cron["DEM_MAX_HS"] == "1.0" and cron["DEM_CELL"] == "2.0" and not notes
              and cfg["STORM_CALM_HS"] == 1.0 and sw.dem_args(cron, True)[:2] == ["--max-hs", "1.0"],
              "waterline_station.conf's DEM_MAX_HS (1.0) wins over the cron's (1.5) for the storm DEMs and the calm "
              "level; its other keys (a forbidden DEM_CELL) are not taken")
        runner = fake_runner_factory()
        ctx = ctx_for(base, NOW, runner=runner)
        ctx.has_hs = True
        E = sw.new_event(ctx, {"events": []}, "obs", "2026-11-14")
        sw.build_dem(ctx, E, sw.evdir(ctx, E) / "pre" / "pre_x", "2026-11-07", "2026-11-13", no_plot=True)
        argv = runner.calls[0] if runner.calls else []
        check("--max-hs" in argv and argv[argv.index("--max-hs") + 1] == "1.0",
              "... and the DEM build is called with --max-hs 1.0")
        (base / "waterline_station.conf").write_text('DEM_MAX_HS="lots"\n')
        cron, notes = sw.cron_settings(base)
        check(cron["DEM_MAX_HS"] == "1.5" and any("not a number" in n for n in notes),
              "a DEM_MAX_HS there that is not a number: the cron's, with a note")

    # 5. recipients: owg.sh's only on CACO-05 in full mode
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "rcpt", settings=None)
        cron, cfg = cfg_for(base)
        to1, w1 = sw.recipients(base, cfg)
        to2, w2 = sw.recipients(base, dict(cfg, STORM_MODE="alerts"))
        (base / "storm_watch_station.conf").write_text(STATION_TEXT)
        to3, w3 = sw.recipients(base, cfg)
        to4, _ = sw.recipients(base, dict(cfg, STORM_EMAIL_TO="x@example.org"))
        check(len(to1) == 3 and w1 == "owg.sh RECIPIENTS" and to2 == [] and "alerts mode" in w2 and to3 == []
              and "storm_watch_station.conf exists" in w3 and to4 == ["x@example.org"],
              "owg.sh's RECIPIENTS only without a station file in full mode; never with a storm_watch_station.conf or "
              "in alerts mode; STORM_EMAIL_TO always")
        start = datetime(2026, 10, 20, tzinfo=UTC)
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_buoy(base / "archive/waves_46026.csv",
                   [(e, h / 0.83 * 1.35, 12.0, 270) for e, h in storm_hours(start, 45, storm_at=storm)])
        sent = []
        days = simple_ground(days_of("2026-10-20", "2026-11-30"))
        run_day(base, "2026-11-12", days)
        with contextlib.redirect_stdout(io.StringIO()):
            ev, _, line = run_day(base, "2026-11-15", days, mailer=lambda *a: sent.append(a) or 0, email=True)
        E = ev.get("2026-11-14_obs") or {}
        check(E and not sent and "WARNING: no email recipients" in line and "no owg.sh fallback" in line
              and [m["status"] for m in E["emails"]] == ["failed"] and evdir_file(base, E, "missed_storm") is not None,
              "a setup station with emails on and no STORM_EMAIL_TO: nothing sent, a WARNING, the body still written")

    # 6. STORM_BUOY / STORM_BUOY2: the fetch, the file names, the GFS-Wave bulletin, the labels
    with tempfile.TemporaryDirectory() as tmp:
        base = station_base(tmp, "buoys", settings='STORM_EMAIL_TO="a@example.org"\n')
        cron, cfg = station_cfg(base)
        net = FakeNet(fixture_routes(), default=(404, b"{}"))
        ctx = sw.Ctx(base, NOW, cfg, cron, mode="live", opener=net, runner=fake_runner_factory(), fetch=True,
                     email=False, local_now=LOCAL, sleep=nosleep)
        state = sw.fresh_state(NOW)
        f = sw.fetch_forecast(ctx, state)
        polar = [u for u, _ in net.seen if "polar.ncep" in u]
        check(polar == ["https://polar.ncep.noaa.gov/waves/WEB/gfswave.latest_run/plots/gfswave.46026.bull"]
              and "gfswave_46026.bull" in f.raw and not any("44008" in k for k in f.raw)
              and sw.source_files(cfg)[-1] == "gfswave_46026.bull",
              "the GFS-Wave bulletin of STORM_BUOY (46026) is fetched and named for it")
        uas = {h.get("User-Agent") for u, h in net.seen if "api.weather.gov" not in u}
        check(uas == {"test01-storm-watch (USGS coastal camera research)"},
              f"the User-Agent names the station (STORM_STATION TEST-01): {uas}")
        with contextlib.redirect_stdout(io.StringIO()):
            rc, line = sw.run_once(ctx)
        arch = sorted(pp.name for pp in (base / "archive/storm_watch/forecasts/2026-11-16").iterdir())
        check(rc == 0 and "1250Z_gfswave_46026.bull.gz" in arch and not any("44008" in n for n in arch),
              f"the archived forecast is named for the buoy ({[n for n in arch if 'gfswave' in n]})")
        raw = {"gfswave_46026.bull": (200, (FIX / "gfswave_44008_SYNTHETIC.bull").read_bytes(), "")}
        fc = sw.build_forecast(raw, NOW, cfg)
        s46 = (fc["groups"].get("gfswave_46026") or [{}])[0]
        check(fc["waves_group"] == "gfswave_46026" and fc["status"].get("gfswave_46026") == "ok"
              and s46.get("name", "").startswith("gfswave 46026 (cycle") and s46.get("cell") == "46026 model point",
              "the bulletin is read under the buoy's name: group, status, series name and cell")
        om = json.loads((FIX / "om_marine_SYNTHETIC.json").read_text())
        om[0]["latitude"], om[0]["longitude"] = 37.75, -122.71
        om[1]["latitude"], om[1]["longitude"] = 37.76, -122.83
        fc = sw.build_forecast({"om_marine.json": (200, json.dumps(om).encode(), "")}, NOW, cfg)
        names = [x["name"] for x in fc["groups"].get("om_marine", [])]
        check(len(names) == 2 and names[0].startswith("om_marine offshore Test (37.75N 122.71W)")
              and names[1].startswith("om_marine 46026 point"), f"the Open-Meteo points labelled for the site and buoy "
                                                               f"({names})")
        cfg_n = sw.load_settings([base / "storm_watch_station.conf"], cron)[0]
        cfg_n2 = dict(cfg_n, STORM_GFSWAVE_STATION=None)
        net2 = FakeNet(fixture_routes(), default=(404, b"{}"))
        ctx2 = sw.Ctx(base, NOW, cfg_n2, cron, mode="dry", opener=net2, fetch=True, email=False, local_now=LOCAL,
                      sleep=nosleep)
        sw.fetch_forecast(ctx2, sw.fresh_state(NOW))
        (base / "storm_watch_settings.conf").write_text("STORM_GFSWAVE_STATION=none\n")
        cron, cfg_none = station_cfg(base)
        check(not [u for u, _ in net2.seen if "polar.ncep" in u] and cfg_none["STORM_GFSWAVE_STATION"] is None
              and sw.gfs_file(cfg_none) is None and len(sw.source_files(cfg_none)) == 6,
              "STORM_GFSWAVE_STATION=none: no bulletin is fetched")
        (base / "storm_watch_settings.conf").write_text("STORM_GFSWAVE_STATION=46013\n")
        cron, cfg_g = station_cfg(base)
        check(sw.gfs_file(cfg_g) == "gfswave_46013.bull" and cfg_g["STORM_BUOY"] == "46026",
              "STORM_GFSWAVE_STATION names another bulletin; the buoy stays STORM_BUOY")
        (base / "storm_watch_settings.conf").write_text("")
        cron, cfg = station_cfg(base)
        # the observed side: archive/waves_46026.csv and waves_46237.csv, never 44008's or 44013's
        t0 = datetime(2026, 11, 10, tzinfo=UTC)
        e = [ep(t0 + timedelta(hours=k)) for k in range(10)]
        write_buoy(base / "archive/waves_46026.csv", [(e[0], 4.0, 12, 270), (e[1], 4.4, 12, 270), (e[2], 4.3, 12, 270)])
        write_buoy(base / "archive/waves_46237.csv", [(e[3], 2.0, 12, 270)])
        write_buoy(base / "archive/waves_44008.csv", [(e[k], 9.0, 12, 270) for k in range(10)])
        write_buoy(base / "archive/waves_44013.csv", [(e[k], 9.0, 12, 270) for k in range(10)])
        o = sw.read_observed(base, t0 + timedelta(hours=9), cfg)
        lab = sw.src_labels(cfg)
        k0 = sw.hour_of(t0) - o.h0
        check(abs(o.prim[k0] - 4.0 * 0.83) < 1e-9 and lab[int(o.psrc[k0])] == "46026 x 0.83"
              and abs(o.prim[k0 + 3] - 2.0 * 0.83) < 1e-9 and lab[int(o.psrc[k0 + 3])] == "46237 converted"
              and np.nanmax(o.raw_hs) == 4.4,
              "the observed waves come from waves_46026.csv (x 0.83) and waves_46237.csv (no transfer file and not "
              "CACO-05: x 0.83, not CACO-05's 44013 fit); 44008's and 44013's files are not read")
        segs = sw.segments(o, cfg)
        check(segs and segs[0]["buoy"] == "46026" and ", 46026 4.4 m" in sw.seg_text(segs[0])
              and "46026 x 0.83" in segs[0]["srcs"], f"a segment names the buoy: {sw.seg_text(segs[0]) if segs else ''}")
        E = {"observed": {"onset": "2026-11-10T00:00:00Z", "peak_hs": 3.65, "peak_time": "2026-11-10T01:00:00Z",
                          "raw_max": 4.4, "end": None, "sources": []}, "storm_dates": None, "forecast": None, "notes": []}
        sl = "\n".join(sw.storm_lines(E, cfg))
        check("46026 up to 4.4 m" in sl and "(Test, nearshore-equivalent)" in sl and "PST" in sl and "44008" not in sl,
              "the email's storm lines: the buoy, the site, local times in STORM_TZ (PST)")
        runner = fake_runner_factory()
        c3 = sw.Ctx(base, NOW, cfg, cron, mode="live", runner=runner, fetch=True, email=False, local_now=LOCAL,
                    sleep=nosleep)
        sw.observe(c3)
        fetched = [(a[a.index("--station") + 1], Path(a[a.index("--output") + 1]).name) for a in runner.calls]
        check(fetched == [("46026", "waves_46026.csv"), ("46237", "waves_46237.csv")],
              f"stale archives: STORM_BUOY and (a setup station, no owg.sh to do it) STORM_BUOY2 fetched ({fetched})")
        net3 = FakeNet([], default=(404, b""))
        with contextlib.redirect_stdout(io.StringIO()):
            rc, _ = sw.climatology(cfg, None, None, "2025-2025", opener=net3)
        check(net3.seen and "/46026h2025.txt.gz" in net3.seen[0][0]
              and net3.seen[0][1]["User-Agent"].startswith("test01-storm-watch"),
              "--climatology defaults to STORM_BUOY's history")
    # --climatology --emit-settings: suggested thresholds as KEY=VALUE lines that --check-settings accepts
    lines = ["#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS  TIDE"]
    t0 = datetime(2025, 1, 1, tzinfo=UTC)
    for k in range(24 * 120):
        tt = t0 + timedelta(hours=k)
        hs = 1.0 + 1.5 * (k % 100) / 99.0 + (2.5 if k % 500 < 6 else 0.0)     # 1.0-2.5 m, a 6 h storm every ~21 days
        lines.append(f"{tt:%Y %m %d %H} 00 MM MM MM {hs:5.2f} 10.00 7.00  60 1010.0 MM MM MM MM MM")
    with tempfile.TemporaryDirectory() as tmp:
        pin = Path(tmp) / "x.txt.gz"
        pin.write_bytes(gzip.compress("\n".join(lines).encode()))
        cfg_f = dict(sw.load_settings(None)[0], STORM_FACING=80.0)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc, res = sw.climatology(cfg_f, [str(pin)], emit=True)
        sug = res.get("suggested") or {}
        emitted = [ln for ln in out.getvalue().splitlines() if re.match(r"^STORM_[A-Z_]+=", ln)]
        base = make_base(tmp, "emit", settings="\n".join(emitted) + "\n")
        rc2, chk = check_settings_run(base)
        check(rc == 0 and sug.get("STORM_DIR_SECTOR") == "0-160" and float(sug["STORM_HS_M"]) >= float(sug["STORM_HS_SUPPORT_M"])
              and float(sug["STORM_CONFIRM_HS"]) < float(sug["STORM_HS_SUPPORT_M"]) and len(emitted) == len(sug) == 8
              and rc2 == 0, f"--climatology --emit-settings prints {len(emitted)} KEY=VALUE lines (STORM_FACING 80 -> "
                            f"STORM_DIR_SECTOR 0-160) that --check-settings accepts ({sug})")
        with contextlib.redirect_stdout(io.StringIO()):
            rc3 = sw.main(["--base", str(base), "--emit-settings"])
        check(rc3 == 1, "--emit-settings without --climatology: refused")
    # the CACO-05 base: a stale 44013 archive is left to owg.sh (no station file)
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "c5fetch")
        runner = fake_runner_factory()
        c4 = ctx_for(base, NOW, fetch=True, runner=runner)
        sw.observe(c4)
        check([a[a.index("--station") + 1] for a in runner.calls] == ["44008"],
              "CACO-05 (no station file): only 44008 is fetched, as today (owg.sh fetches 44013)")

    # 7. the Open-Meteo sea-cell tests: STORM_FACING's sector; an empty STORM_MIN_CELL_LON turns the old test off
    om = json.loads((FIX / "om_marine_SYNTHETIC.json").read_text())
    pts = [(41.90, -69.80), (40.50, -69.25)]
    bay = json.loads(json.dumps(om))
    bay[0]["longitude"] = -70.05
    r_old = sw.parse_om_marine(bay, [(41.90, -69.90), pts[1]], NOW, 96, -69.93)
    r_off = sw.parse_om_marine(bay, [(41.90, -69.90), pts[1]], NOW, 96, None)
    r_fac = sw.parse_om_marine(bay, [(41.90, -69.90), pts[1]], NOW, 96, None, facing=80.0, beach=(41.89, -69.96))
    r_ok = sw.parse_om_marine(om, pts, NOW, 96, None, facing=80.0, beach=(41.89, -69.96))
    check(isinstance(r_old[0], str) and isinstance(r_off[0], dict) and isinstance(r_fac[0], str)
          and "outside the sector it faces (80 +- 80 deg)" in r_fac[0] and isinstance(r_ok[0], dict)
          and isinstance(r_fac[1], dict),
          "a cell in Cape Cod Bay: refused by the longitude test; accepted with STORM_MIN_CELL_LON empty; refused by "
          "STORM_FACING 80 (it lies at ~270 deg from the beach); the offshore cell (85 deg) accepted; the buoy point "
          "is never tested")
    west = json.loads(json.dumps(om))
    west[0]["latitude"], west[0]["longitude"] = 37.75, -122.70
    behind = json.loads(json.dumps(om))
    behind[0]["latitude"], behind[0]["longitude"] = 37.75, -122.40
    pw = [(37.75, -122.60), pts[1]]
    rw = sw.parse_om_marine(west, pw, NOW, 96, None, facing=270.0, beach=(37.75, -122.51))
    rb = sw.parse_om_marine(behind, pw, NOW, 96, None, facing=270.0, beach=(37.75, -122.51))
    check(isinstance(rw[0], dict) and isinstance(rb[0], str) and "outside the sector" in rb[0],
          "a west-facing beach (270): a cell to the west kept, one behind the beach (east) refused")

    # 8. STORM_TZ: zoneinfo first, else pandas (Python 3.8, or a zoneinfo that cannot be imported)
    t = datetime(2026, 11, 17, 14, tzinfo=UTC)
    try:
        import zoneinfo     # noqa: F401
        have_zi = True
    except ImportError:
        have_zi = False
    z = sw.to_zone(t, "America/Los_Angeles")
    a_txt = sw.local_time(t, "America/Los_Angeles")
    saved = sys.modules.get("zoneinfo")
    had = "zoneinfo" in sys.modules
    sys.modules["zoneinfo"] = None      # 'import zoneinfo' now raises ImportError, as on Python 3.8
    try:
        z2 = sw.to_zone(t, "America/Los_Angeles")
        b_txt = sw.local_time(t, "America/Los_Angeles")
        valid = (sw.valid_tz("America/Los_Angeles"), sw.valid_tz("Mars/Olympus_Mons"), sw.valid_tz("../etc"))
        e_txt = sw.when(t, "Pacific/Honolulu")
    finally:
        if had:
            sys.modules["zoneinfo"] = saved
        else:
            del sys.modules["zoneinfo"]
    check((type(z.tzinfo).__name__ == "ZoneInfo") == have_zi and a_txt == "Tue 17 Nov 06:00 PST",
          f"STORM_TZ through zoneinfo when Python has it ({'yes' if have_zi else 'no: pandas'}): {a_txt}")
    check(type(z2).__name__ == "Timestamp" and b_txt == a_txt and valid == (True, False, False)
          and e_txt == "17 Nov 14Z (Tue 17 Nov 04:00 HST)",
          f"... and through pandas when zoneinfo cannot be imported: {b_txt}; {e_txt}; zone check {valid}")
    check(sw.local_time(t) == "Tue 17 Nov 09:00 EST" and sw.when(t) == "17 Nov 14Z (Tue 17 Nov 09:00 EST)",
          "no STORM_TZ: CACO-05's America/New_York, as today")
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "s.conf"
        p.write_text("STORM_TZ=Mars/Olympus_Mons\nSTORM_BUOY=../x\nSTORM_MODE=sometimes\n")
        c, notes = sw.load_settings(p)
        check(c["STORM_TZ"] == "America/New_York" and c["STORM_BUOY"] == "44008" and c["STORM_MODE"] == "full"
              and len(notes) == 3, "a bad zone, buoy id or mode in a running station's file: the default, with a note")


# --------------------------------------------------------------------------------------------
# T15 alerts mode (STORM_MODE=alerts): no DEM; a storm-over email; closed 72 h after the end; events.csv
# --------------------------------------------------------------------------------------------

def events_rows(base):
    p = base / "archive/storm_watch/events.csv"
    if not p.exists():
        return []
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def bodies(base, E):
    return {m["kind"]: (base / "archive/storm_watch/events" / E["id"] / m["body_file"]).read_text()
            for m in E["emails"] if m.get("body_file")}


def t15():
    print("T15 alerts mode: no DEM, storm over, closed 72 h after the end, events.csv")
    start = datetime(2026, 10, 20, tzinfo=UTC)
    gi = simple_ground(days_of("2026-10-20", "2026-11-30"))
    dem_word = re.compile(r"\bDEMs?\b|NAVD88")
    # 1. a storm the forecast missed, on CACO-05's defaults in alerts mode
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "alerts", settings='STORM_MODE=alerts\nSTORM_EMAIL_TO="a@example.org"\n')
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(start, 45, storm_at=storm)])
        runner = fake_runner_factory()
        sent = []
        got = {}
        for d in days_of("2026-11-12", "2026-11-20"):
            ev, ctx, line = run_day(base, d, gi, runner=runner, mailer=lambda s, t, b, a: sent.append((s, t, b)) or 0,
                                    email=True)
            E = ev.get("2026-11-14_obs")
            got[d] = (E["status"], E["closed_reason"], [m["kind"] for m in E["emails"]]) if E else None
        E = ev["2026-11-14_obs"]
        end = sw.parse_iso(E["observed"]["end"])
        check(got["2026-11-15"] == ("storm", None, ["missed_storm"])
              and got["2026-11-16"] == ("post", None, ["missed_storm", "storm_over"])
              and got["2026-11-17"][0] == "post" and got["2026-11-18"][:2] == ("closed", "complete")
              and end == datetime(2026, 11, 15, 6, tzinfo=UTC),
              f"the storm ends 15 Nov 06Z: a storm-over email on the run that sees the end (16 Nov), still open on "
              f"17 Nov, closed on the first run 72 h after the end (18 Nov) ({got})")
        check(any("72 h after the storm's end" in h for h in E["history"]) and E["pre"]["status"] == "planned"
              and E["post3"] is None and E["post7"] is None and not E["change"] and E.get("mode") == "alerts",
              "... closed by the clock: no pre, post or change step ran")
        check(not any(Path(a[1]).name in ("dem_from_contours.py", "dem_change.py", "dem_map.py", "asc_to_geotiff.py")
                      for a in runner.calls)
              and not (base / "archive/storm_watch/work").exists()
              and json.loads((base / "archive/storm_watch/state.json").read_text()).get("ground_cache") is None,
              "no DEM, map or change call; the ground file is never read (no working copy, no ground cache)")
        rows = events_rows(base)
        check(len(rows) == 1 and rows[0]["event"] == "2026-11-14_obs" and rows[0]["closed_reason"] == "complete"
              and rows[0]["storm_end_utc"] == "2026-11-15T06:00:00Z" and rows[0]["mode"] == "alerts"
              and rows[0]["closed_utc"].startswith("2026-11-18") and rows[0]["buoy"] == "44008"
              and rows[0]["emails"] == "missed_storm:sent storm_over:sent" and "@" not in json.dumps(rows),
              f"one events.csv row, written by the closing run, never again ({rows[:1]})")
        texts = bodies(base, E)
        check(set(texts) == {"missed_storm", "storm_over"} and not any(dem_word.search(t) for t in texts.values())
              and all("alerts mode" in t.lower() for t in texts.values())
              and [s[:42] for s, _, _ in sent] == ["[CACO-05 storm watch] A storm the forecast",
                                                     "[CACO-05 storm watch] Storm over at Marcon"],
              "the emails (missed storm, storm over) carry no DEM or NAVD88 wording; both sent")
        check("is over: the waves at Marconi stayed below 1.5 m for 24 hours" in texts["storm_over"]
              and "closes 72 h after the end" in texts["storm_over"], "the storm-over email says why and when it closes")
    # 2. a trigger for the same storm while the post event is open (72 h) belongs to it: no second heads-up
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "tail", settings='STORM_MODE=alerts\nSTORM_EMAIL=0\n')
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(start, 45, storm_at=storm)])
        for d in ("2026-11-12", "2026-11-15", "2026-11-16"):
            run_day(base, d, gi)

        def stale_fc(c):
            # an old model run still carrying the storm: rough from 15 Nov 00Z (before its end) to 18 Nov
            h0 = sw.hour_of(datetime(2026, 11, 15, 0, tzinfo=UTC))
            hs = np.where(np.arange(120) < 66, 3.8, 1.0)
            return fc_of([series(hs, h0=h0)])
        ev, _, _ = run_day(base, "2026-11-17", gi, stale_fc)
        check(list(ev) == ["2026-11-14_obs"] and any("forecast still triggering" in n for n in ev["2026-11-14_obs"]["notes"]),
              "a forecast for the storm's tail within the 72 h: noted on the open event, no new event or heads-up")
    # 3. a forecast storm that never comes: a false alarm without a pending pre DEM; a setup station's names and zone
    with tempfile.TemporaryDirectory() as tmp:
        base = station_base(tmp, "fa", extra="STORM_MODE=alerts\nSTORM_EMAIL=0\n", settings="")
        write_buoy(base / "archive/waves_46026.csv",
                   [(e, 1.2, 10.0, 270) for e, _ in storm_hours(start, 45)])
        cron, cfg = station_cfg(base)

        def west_storm(c):
            fc = storm_fc(c, datetime(2026, 11, 13, 18, tzinfo=UTC))
            fc["waves"][0]["dir"][:] = 270.0
            return fc
        runner = fake_runner_factory()
        for d in days_of("2026-11-10", "2026-11-18"):
            ev, ctx, line = run_day(base, d, gi, west_storm if d == "2026-11-10" else None, runner=runner)
        E = list(ev.values())[0]
        texts = bodies(base, E)
        check(E["closed_reason"] == "false_alarm" and not E["pending"] and set(texts) == {"heads_up", "false_alarm"}
              and not any(dem_word.search(t) for t in texts.values()),
              f"alerts mode, a forecast storm that never comes: a false alarm, no pre DEM left pending, no DEM wording "
              f"({E['closed_reason']}, {list(texts)})")
        check("[TEST-01 storm watch] Storm forecast for Test Beach: onset ~" in (
              base / "archive/storm_watch/events" / E["id"] / E["emails"][0]["body_file"]).read_text()
              and "PST" in texts["heads_up"] and "alerts mode" in texts["heads_up"]
              and "(46026, 46237)" in texts["heads_up"],
              "the heads-up: the station's subject and site, local times in its zone, its buoys")
        rows = events_rows(base)
        check(len(rows) == 1 and rows[0]["closed_reason"] == "false_alarm" and rows[0]["buoy"] == "46026"
              and rows[0]["station"] == "TEST-01" and rows[0]["forecast_rules"] == "a"
              and not any(Path(a[1]).name == "dem_from_contours.py" for a in runner.calls),
              "its events.csv row; no DEM built")
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            sw.main(["--base", str(base), "--status"])
        check("alerts mode" in out.getvalue() and "     pre " not in out.getvalue(),
              "--status says alerts mode, no pre line")
    # 4. full mode (the default) writes no events.csv and keeps its emails as today
    with tempfile.TemporaryDirectory() as tmp:
        base = make_base(tmp, "full")
        storm = datetime(2026, 11, 14, 3, tzinfo=UTC)
        write_marconi(base / "archive/waves_marconi.csv",
                      [(e, h, "blend", 1, 60, 9) for e, h in storm_hours(start, 45, storm_at=storm)])
        for d in ("2026-11-12", "2026-11-15", "2026-11-16", "2026-11-19"):
            ev, _, line = run_day(base, d, gi)
        E = ev["2026-11-14_obs"]
        check(not events_rows(base) and "storm_over" not in [m["kind"] for m in E["emails"]] and "mode" not in E
              and "csv_row" not in E and "alerts mode" not in line and E["status"] == "post",
              "full mode: no storm-over email, no events.csv, no new keys in the event, nothing closed by the clock")


# --------------------------------------------------------------------------------------------
# T16 storm_watch.sh on a setup-managed NUC; the station file's exports
# --------------------------------------------------------------------------------------------

def t16():
    print("T16 storm_watch.sh: a setup-managed NUC's crontab; waterline_station.conf's exports")
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp) / "station"
        base.mkdir()
        shutil.copy(str(HERE / "storm_watch.sh"), str(base))
        b = Path(tmp) / "bin"
        b.mkdir()
        tab = Path(tmp) / "crontab.txt"
        tab.write_text("30 12 * * * /mnt/I2Rgus_Data/waterline/waterline_timex_cron.sh\n")
        (b / "crontab").write_text(f'#!/bin/bash\nif [ "$1" = "-l" ]; then cat "{tab}"; else cat > "{tab}"; fi\n')
        (b / "crontab").chmod(0o755)
        home = Path(tmp) / "home"
        (home / "stationTools").mkdir(parents=True)
        (home / "stationTools/crontab_new_station.txt").write_text("# the crontab setup installed\n")
        env = dict(os.environ, PATH=f"{b}:{os.environ['PATH']}", HOME=str(home))
        sh = str(base / "storm_watch.sh")
        before = tab.read_text()
        r1 = subprocess.run([sh, "--install-cron"], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
        r2 = subprocess.run([sh, "--remove-cron"], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
        check(r1.returncode == 1 and 'STORM_WATCH="yes"' in r1.stdout and "./setup.sh --only 10" in r1.stdout
              and r2.returncode == 1 and 'STORM_WATCH="no"' in r2.stdout and tab.read_text() == before
              and not (base / "logs/crontab.before_storm_watch_install").exists(),
              "~/stationTools/crontab_new_station.txt exists: --install-cron and --remove-cron refuse, point to "
              "STORM_WATCH in station.conf and ./setup.sh --only 10, and leave the crontab as it was")
        env2 = dict(env, HOME=str(Path(tmp) / "elsewhere"))
        r3 = subprocess.run([sh, "--install-cron"], env=env2, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True)
        check(r3.returncode == 0 and "storm_watch.sh" in tab.read_text(), "... and without it, --install-cron works")
        # the station file: checked (a fake waterline_station_check.py exit), then sourced; its exports reach python
        envout = Path(tmp) / "env.txt"
        (b / "python3").write_text('#!/bin/bash\ncase "$1" in\n  *waterline_station_check.py) echo "check $2"; '
                                   'exit "${CHECK_RC:-0}" ;;\nesac\nenv > "$ENV_OUT"\nexit 0\n')
        (b / "python3").chmod(0o755)
        (base / "waterline_station.conf").write_text('STATION_NAME="TEST-01"\nGEOTIFF_EPSG="26910"\nCAMS="c1"\n'
                                                     'export STATION_NAME GEOTIFF_EPSG\n')
        env3 = dict(env2, ENV_OUT=str(envout))
        env3.pop("STATION_NAME", None)
        env3.pop("GEOTIFF_EPSG", None)
        r4 = subprocess.run([sh], env=env3, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                            universal_newlines=True)
        seen = envout.read_text().splitlines() if envout.exists() else []
        check(r4.returncode == 0 and "STATION_NAME=TEST-01" in seen and "GEOTIFF_EPSG=26910" in seen
              and not any(ln.startswith("CAMS=") for ln in seen)
              and "check check-overlay" in (base / "logs/storm_watch_detail.log").read_text(),
              "waterline_station.conf: check-overlay first, then sourced; its exported keys reach storm_watch.py "
              "(STATION_NAME, GEOTIFF_EPSG), the others do not")
        envout.unlink()
        r5 = subprocess.run([sh], env=dict(env3, CHECK_RC="1"), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            stdin=subprocess.DEVNULL, universal_newlines=True)
        runlog = (base / "logs/storm_watch.log").read_text()
        check(r5.returncode == 1 and not envout.exists() and "ERROR: waterline_station.conf refused" in runlog,
              "a station file check-overlay refuses: storm_watch.py is not run, an ERROR line, exit 1")
        r6 = subprocess.run([sh, "--dry-run"], env=dict(env3, CHECK_RC="1"), stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, universal_newlines=True)
        check(r6.returncode == 0 and envout.exists() and "STATION_NAME=TEST-01" not in envout.read_text(),
              "--dry-run (builds nothing) neither checks nor sources it")


# --------------------------------------------------------------------------------------------
# real captures (on the NUC)
# --------------------------------------------------------------------------------------------

def real_captures(folder):
    print(f"REAL captures in {folder} ([NUC])")
    cfg = sw.load_settings(None)[0]
    d = Path(folder)
    raw, stamps = {}, []
    for name in sw.SOURCE_FILES:
        p = d / name
        if p.exists():
            raw[name] = (200, p.read_bytes(), "")
            stamps.append(p.stat().st_mtime)
    now = datetime.fromtimestamp(max(stamps), tz=UTC) if stamps else datetime.now(UTC)   # when they were saved
    if not raw:
        check(False, "no capture found (run ./storm_watch.sh --dry-run --save-raw DIR first)")
        return
    fc = sw.build_forecast(raw, now, cfg)
    for k, v in sorted(fc["status"].items()):
        check(v.startswith(("ok", "no wind")) or "wind only" in v, f"{k}: {v}")
    if "nws_points.json" in raw:
        p = sw.parse_nws_points(json.loads(raw["nws_points.json"][1]))
        check(p["forecastGridData"].startswith("https://api.weather.gov/gridpoints/"),
              f"/points: {p['gridId']} {p['gridX']},{p['gridY']}, zone {str(p['forecastZone']).rsplit('/', 1)[-1]}")
    for s in fc["waves"] + fc["waves_other"]:
        check(np.nanmax(s["hs"]) < 20 and np.isfinite(s["hs"]).sum() >= 24, f"{s['name']}: Hs up to {np.nanmax(s['hs']):.1f} m")
    if fc["alerts"] is not None:
        print(f"  alerts kept: {[(a['code'], a['event'], a['zones']) for a in fc['alerts']]}")
    T = sw.evaluate_trigger(fc, now, cfg)
    print(f"  trigger: {T['triggered']} {T['rules']}; peak {T['peak']}")


def main():
    if "--real-captures" in sys.argv:
        real_captures(sys.argv[sys.argv.index("--real-captures") + 1])
    else:
        os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="mpl_storm_watch_"))
        for t in (t1, t2, t3, t4, t5, t6, t7, t8, t9, t9b, t10, t11, t12, t13, t14, t15, t16):
            t()
    if FAILED:
        print(f"{len(FAILED)} storm watch check(s) FAILED:")
        for f in FAILED:
            print("  " + f)
        raise SystemExit(1)
    print("all storm watch tests passed")


if __name__ == "__main__":
    main()

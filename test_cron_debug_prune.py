#!/usr/bin/env python3
"""
Tests For The Station Cron's Debug-Image Pruning
===================================================
Plain python3 (no pytest needed on the station):

    python3 test_cron_debug_prune.py [--keep DIR]

WHAT. waterline_timex_cron.sh removes the detector's debug images (one
folder per frame, ~8 MB each) DEBUG_KEEP_DAYS days after they were written
(prune_debug_images). An rm -rf on a configurable folder must remove what
it is meant to and nothing else:
  * the detector's frame folders (named after their photo: a 10-digit
    Argus epoch and a dot first) older than the limit go, names with
    spaces too; younger ones stay;
  * anything else in the folder stays: a folder kept there by hand, a
    plain file, a link to a folder elsewhere (and what it points to), and
    what a link inside a removed frame folder points to;
  * DEBUG_DIR may itself be a link to another disk: its old frame folders
    go (before, find looked at the link and pruned nothing);
  * DEBUG_KEEP_DAYS="" keeps everything and says nothing; a value that is
    not a whole number is a WARNING in the log and removes nothing; 08 is
    8 days (not a bad octal number);
  * a folder that could not be removed is counted as left, with a WARNING,
    not as removed -- also one rm emptied but could not remove (emptying
    it made it new again, and a second find no longer counted it);
  * a folder reaching the limit while rm runs is neither removed nor
    counted (one cut-off time for every find): no false WARNING.

HOW. The function is taken from the cron's own text (so the test runs what
the station runs) and called by bash on a scratch tree in a temporary
folder (--keep DIR to look at it), with log() writing to a file.

STATION FILES. None are read.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
CRON = HERE / "waterline_timex_cron.sh"

FAILED = []


def check(cond, msg):
    print(("  ok  " if cond else "  FAIL ") + msg)
    if not cond:
        FAILED.append(msg)


def prune_function():
    """The text of prune_debug_images() in the cron, from its first line to the closing brace."""
    m = re.search(r"^prune_debug_images\(\) \{\n.*?^\}\n", CRON.read_text(), re.S | re.M)
    return m.group(0) if m else ""


def run_prune(d, folder, keep, path_prefix=None):
    """Calls the cron's prune_debug_images on `folder` with `keep` days. -> (exit status, log text)."""
    logf = Path(d) / "prune.log"
    if logf.exists():
        logf.unlink()
    script = Path(d) / "prune_driver.sh"
    script.write_text("set -u\n"
                      'log() { echo "$*" >> "$LOGF"; }\n'
                      + prune_function()
                      + 'prune_debug_images "$1" "$2"\n')
    env = dict(os.environ, LOGF=str(logf))
    if path_prefix:
        env["PATH"] = f"{path_prefix}:{env['PATH']}"
    r = subprocess.run(["bash", str(script), str(folder), keep], env=env, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True)
    return r.returncode, (logf.read_text() if logf.exists() else "") + r.stdout


def age(path, days):
    """Sets the modification time of `path` (not following a link) to `days` days ago."""
    t = time.time() - days * 86400
    os.utime(str(path), (t, t), follow_symlinks=False)


FRAME_OLD = "1727712000.Mon.Sep.30_16_00_00.GMT.2024.CACO05.c1.timex"
FRAME_SPACES = "1727713800.Mon Sep 30 16_30_00.GMT.2024.CACO05.c2.timex"
FRAME_LINKED = "1727715600.Mon.Sep.30_17_00_00.GMT.2024.CACO05.c1.timex"
FRAME_NEW = "1728316800.Mon.Oct.07_16_00_00.GMT.2024.CACO05.c1.timex"


def tree(d):
    """A debug folder as the detector leaves it, plus what must survive. -> (debug folder, outside folder)."""
    dbg = Path(d) / "debug"
    out = Path(d) / "outside"
    (out / "keepme").mkdir(parents=True)
    (out / "keepme" / "precious.txt").write_text("must survive\n")
    for name in (FRAME_OLD, FRAME_SPACES, FRAME_LINKED, FRAME_NEW, "kept by hand"):
        (dbg / name).mkdir(parents=True)
        (dbg / name / "energy_final.png").write_bytes(b"png")
    (dbg / FRAME_LINKED / "to_outside").symlink_to(out / "keepme")            # a link inside a frame folder
    (dbg / "1727717400.Mon.Sep.30_17_30_00.GMT.2024.CACO05.c2.timex").symlink_to(out / "keepme")  # a link
    (dbg / "1727719200.loose.png").write_bytes(b"png")                         # a plain file
    for name in (FRAME_OLD, FRAME_SPACES, FRAME_LINKED, "kept by hand", "1727719200.loose.png",
                 "1727717400.Mon.Sep.30_17_30_00.GMT.2024.CACO05.c2.timex"):
        age(dbg / name, 10)
    age(dbg / FRAME_NEW, 6)
    age(out / "keepme", 10)
    return dbg, out


def test_prune(d):
    print("prune_debug_images: what goes and what stays")
    d = Path(d) / "plain"
    d.mkdir()
    dbg, out = tree(d)
    rc, log = run_prune(d, dbg, "7")
    left = sorted(p.name for p in dbg.iterdir())
    check(rc == 0 and left == sorted([FRAME_NEW, "kept by hand", "1727719200.loose.png",
                                      "1727717400.Mon.Sep.30_17_30_00.GMT.2024.CACO05.c2.timex"]),
          f"the 3 frame folders 10 days old went (one with spaces, one holding a link); left: {left}")
    check((out / "keepme" / "precious.txt").read_text() == "must survive\n",
          "what the links pointed to is untouched")
    check("removed 3 frame folder(s) 7 or more days old" in log and "WARNING" not in log,
          f"logged as removed, no WARNING ({log.strip()!r})")
    rc, log = run_prune(d, dbg, "7")
    check(rc == 0 and log == "", "a second run: nothing old left, nothing logged")


def test_linked_folder(d):
    print("prune_debug_images: DEBUG_DIR a link to another disk")
    d = Path(d) / "linked"
    d.mkdir()
    dbg, out = tree(d / "disk")
    (d / "base").mkdir()
    (d / "base" / "debug").symlink_to(dbg)
    rc, log = run_prune(d, d / "base" / "debug", "7")
    left = sorted(p.name for p in dbg.iterdir())
    check(rc == 0 and FRAME_OLD not in left and FRAME_NEW in left and "kept by hand" in left
          and "removed 3 frame folder(s)" in log and (d / "base" / "debug").is_symlink(),
          f"the old frame folders behind the link went, the link and the rest stayed ({log.strip()!r})")


def test_settings(d):
    print("prune_debug_images: DEBUG_KEEP_DAYS empty, or not a number")
    d = Path(d) / "settings"
    d.mkdir()
    dbg, _ = tree(d)
    before = sorted(p.name for p in dbg.iterdir())
    rc, log = run_prune(d, dbg, "")
    check(rc == 0 and log == "" and sorted(p.name for p in dbg.iterdir()) == before,
          "DEBUG_KEEP_DAYS='': everything kept, nothing logged")
    for bad in ("7d", "-3", "seven"):
        rc, log = run_prune(d, dbg, bad)
        check(rc == 0 and log.startswith(f"WARNING: DEBUG_KEEP_DAYS='{bad}' is not a whole number")
              and sorted(p.name for p in dbg.iterdir()) == before,
              f"DEBUG_KEEP_DAYS='{bad}': a WARNING, nothing removed, the job carries on")
    rc, log = run_prune(d, d / "missing", "7")
    check(rc == 0 and log == "", "no debug folder yet: nothing to do")
    # a leading zero is still days, in base 10 (bash's arithmetic read 08 as a bad octal number, and
    # the function stopped there: nothing pruned, nothing logged)
    rc, log = run_prune(d, dbg, "08")
    left = sorted(p.name for p in dbg.iterdir())
    check(rc == 0 and FRAME_OLD not in left and FRAME_NEW in left and "kept by hand" in left
          and "removed 3 frame folder(s) 08 or more days old" in log and "WARNING" not in log,
          f"DEBUG_KEEP_DAYS='08': 8 days, the 10-day-old frame folders went ({log.strip()!r})")


def test_not_removed(d):
    print("prune_debug_images: folders that could not be removed")
    d = Path(d) / "stuck"
    d.mkdir()
    dbg, _ = tree(d)
    fake = d / "bin"
    fake.mkdir()
    (fake / "rm").write_text("#!/bin/sh\nexit 1\n")       # an rm that removes nothing (as without permission)
    (fake / "rm").chmod(0o755)
    rc, log = run_prune(d, dbg, "7", path_prefix=str(fake))
    check(rc == 0 and FRAME_OLD in {p.name for p in dbg.iterdir()}
          and "removed 0 frame folder(s)" in log
          and "WARNING: debug images: 3 frame folder(s) 7 or more days old could not be removed" in log,
          f"counted again afterwards: 0 removed, a WARNING for the 3 left ({log.strip()!r})")
    # an rm that empties each folder but cannot remove the folder itself (DEBUG_DIR not writable by
    # the job, or a subfolder of another user): emptying a folder makes it new again, and a second
    # find counting the old folders no longer saw it -- logged as removed, no WARNING
    d = d.parent / "emptied"
    d.mkdir()
    dbg, out = tree(d)
    fake = d / "bin"
    fake.mkdir()
    (fake / "rm").write_text('#!/bin/sh\nshift\nfor p in "$@"; do find "$p" -mindepth 1 -delete; done\nexit 1\n')
    (fake / "rm").chmod(0o755)
    rc, log = run_prune(d, dbg, "7", path_prefix=str(fake))
    left = {p.name for p in dbg.iterdir()}
    check(rc == 0 and {FRAME_OLD, FRAME_SPACES, FRAME_LINKED} <= left
          and not any((dbg / n / "energy_final.png").exists() for n in (FRAME_OLD, FRAME_SPACES))
          and (out / "keepme" / "precious.txt").read_text() == "must survive\n"
          and "removed 0 frame folder(s)" in log
          and "WARNING: debug images: 3 frame folder(s) 7 or more days old could not be removed" in log,
          f"emptied but still there: 0 removed, a WARNING for the 3 left ({log.strip()!r})")


def test_limit_reached_during_rm(d):
    print("prune_debug_images: a folder reaching the limit while rm runs")
    d = Path(d) / "race"
    d.mkdir()
    dbg = d / "debug"
    edge = "1728316800.Mon.Oct.07_16_00_00.GMT.2024.CACO05.c2.timex"
    for name in (FRAME_OLD, edge):
        (dbg / name).mkdir(parents=True)
    fake = d / "bin"
    fake.mkdir()
    (fake / "rm").write_text(f'#!/bin/sh\nsleep 3\nexec {shutil.which("rm")} "$@"\n')   # a slow rm
    (fake / "rm").chmod(0o755)
    age(dbg / FRAME_OLD, 10)
    age(dbg / edge, 7 - 2 / 86400)          # 2 s short of 7 days when the function starts, 7 days after rm
    rc, log = run_prune(d, dbg, "7", path_prefix=str(fake))
    left = sorted(p.name for p in dbg.iterdir())
    # each find measured the age from its own start, so a recount after the 3 s rm found the edge
    # folder old, never tried: 'removed 0' and a WARNING (the station status WARN for that run)
    check(rc == 0 and left == [edge] and "removed 1 frame folder(s) 7 or more days old" in log
          and "WARNING" not in log,
          f"one cut-off for every find: the old folder removed, the other left for the next run, "
          f"no WARNING ({log.strip()!r})")


def test_cron_calls_it():
    print("waterline_timex_cron.sh: defines it before step 1b calls it, with its settings")
    text = CRON.read_text()
    f, c = text.find("prune_debug_images() {"), text.find('prune_debug_images "$DEBUG_DIR" "$DEBUG_KEEP_DAYS"')
    check(prune_function() != "" and 0 <= f < c, "the function, then the call after the detector")
    check(text.find('--debug-dir "$DEBUG_DIR"') < c, "the call comes after the detector has written this run's images")
    r = subprocess.run(["bash", "-n", str(CRON)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True)
    check(r.returncode == 0, f"bash -n {CRON.name} {r.stdout.strip()}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--keep", default=None, help="keep the outputs, in a new folder inside this one")
    args = ap.parse_args()
    if args.keep:
        Path(args.keep).mkdir(parents=True, exist_ok=True)
    d = os.path.realpath(tempfile.mkdtemp(prefix="cron_prune_", dir=args.keep))
    print(f"test folder: {d}")
    try:
        test_cron_calls_it()
        test_prune(d)
        test_linked_folder(d)
        test_settings(d)
        test_not_removed(d)
        test_limit_reached_during_rm(d)
    finally:
        if not args.keep:
            shutil.rmtree(d, ignore_errors=True)
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED")
        return 1
    print("all cron debug-prune tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

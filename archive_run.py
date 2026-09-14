#!/usr/bin/env python3
"""
Archive Run
-------------
Copies the CURRENT contents of ground_truth/ and processed/ into a
permanent, dated archive folder that survives both the nightly source
image wipe and each detector run's own cleanup (waterline_detector_v5.py
clears input/processed at the start of every run).

WHY THIS EXISTS: every aggregation tool built for this project
(derive_envelope.py, derive_bias_correction.py, batch_bias_analysis.py)
can only see whatever happens to still be sitting in ground_truth/ and
processed/ RIGHT NOW. Source images are cleared nightly, and
processed/ gets wiped at the start of every detector run -- so a full
day's ground-truth-to-detector-output pairing, once made, has no home
that survives past the next run. This is why "how much data do we
actually have" kept changing between sessions: July 24's ground truth
was correctly paired and analyzed multiple times, but the matching
processed/ output was never preserved anywhere, so by the time a
broader multi-day correction was needed, that day's detector output
was already gone for good (the source images that would let it be
regenerated are wiped nightly too).

Run this AFTER each detector run + ground truth collection, BEFORE the
next run's cleanup. It copies (not moves) so today's ground_truth/ and
processed/ stay usable for immediate follow-up analysis, while also
building up a permanent, ever-growing record that derive_envelope.py /
derive_bias_correction.py / batch_bias_analysis.py can eventually be
pointed at for a real multi-day, multi-season aggregate.

Usage:
    python3 archive_run.py [optional label, e.g. "2026-07-27-morning"]
"""

import sys
import shutil
from pathlib import Path
from datetime import datetime, timezone


GROUND_TRUTH_FOLDER = "/mnt/I2Rgus_Data/waterline/ground_truth"
PROCESSED_FOLDER = "/mnt/I2Rgus_Data/waterline/processed"
ARCHIVE_FOLDER = "/mnt/I2Rgus_Data/waterline/archive"


def copy_new_files(src_dir, dst_dir):
    """
    Copies every file from src_dir into dst_dir, skipping files that
    already exist there with the same name (so re-running this
    multiple times, or archiving overlapping batches, doesn't
    error out or churn -- it just fills in whatever's new).
    """
    src_dir = Path(src_dir)
    dst_dir = Path(dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    copied, skipped = 0, 0
    for f in src_dir.glob("*"):
        if not f.is_file():
            continue
        target = dst_dir / f.name
        if target.exists():
            skipped += 1
            continue
        shutil.copy2(f, target)
        copied += 1
    return copied, skipped


def main():
    label = sys.argv[1] if len(sys.argv) > 1 else datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")

    archive_gt = Path(ARCHIVE_FOLDER) / "ground_truth"
    archive_processed = Path(ARCHIVE_FOLDER) / "processed"

    print(f"Archiving current ground_truth/ and processed/ (label: {label})")
    print(f"  Source ground truth : {GROUND_TRUTH_FOLDER}")
    print(f"  Source processed    : {PROCESSED_FOLDER}")
    print(f"  Archive destination : {ARCHIVE_FOLDER}")
    print()

    gt_copied, gt_skipped = copy_new_files(GROUND_TRUTH_FOLDER, archive_gt)
    print(f"Ground truth : {gt_copied} new file(s) archived, {gt_skipped} already present (skipped)")

    proc_copied, proc_skipped = copy_new_files(PROCESSED_FOLDER, archive_processed)
    print(f"Processed    : {proc_copied} new file(s) archived, {proc_skipped} already present (skipped)")

    print()
    print("Done. This archive is permanent and NOT touched by the detector's")
    print("own cleanup or the nightly source image wipe. Point derive_envelope.py /")
    print("derive_bias_correction.py / batch_bias_analysis.py at the archive")
    print("folders (instead of the working ground_truth/processed folders)")
    print("once you want a real multi-day/multi-season aggregate -- the working")
    print("folders will keep reflecting only whatever the most recent run left behind.")


if __name__ == "__main__":
    main()

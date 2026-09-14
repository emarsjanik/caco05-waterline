#!/usr/bin/env python3
"""
Batch Ground Truth Collection Tool (v3)
------------------------------------------
Loops through images in the detector's input folder, opening each one
for you to click the true shoreline.

WHAT CHANGED IN v3 (and why):
  * MOUSE-ONLY CONTROLS. v2 relied on keyboard shortcuts (Enter / s /
    q). On this station's matplotlib + GTK3 backend, key events whose
    keysym can't be mapped raise inside matplotlib's OWN handler
    ("AttributeError: 'NoneType' object has no attribute 'lower'"),
    which aborts the callback chain BEFORE the script's key handler
    runs -- so Enter appeared to do nothing while printing tracebacks.
    There are now on-screen buttons (Save & Next / Skip / Undo /
    Quit), so collection never depends on the keyboard. Keyboard
    shortcuts still work where the backend allows, but are no longer
    required.
  * KEYSYM CRASH GUARDED. The underlying matplotlib bug is patched
    defensively at import (see _patch_gtk3_keysym_bug) so the
    tracebacks stop cluttering the console even if you do press keys.
  * RESUME. Images that already have a ground-truth CSV are skipped
    automatically, so a 58-image session can be done across several
    sittings without redoing work or accidentally overwriting it.
    Use --redo to override.
  * FILTERING. --filter SUBSTRING restricts collection to matching
    filenames, e.g. --filter 1789 for one day's captures, so you can
    prioritise current imagery instead of working through everything.

Controls per image:
    Left click    : add a point
    Right click   : undo last point
    [Save & Next] : save this image's line and advance
    [Skip]        : no ground truth saved for this image
    [Undo]        : remove last point
    [Quit]        : stop, keeping everything collected so far

Usage:
    python3 collect_all_ground_truth.py
    python3 collect_all_ground_truth.py --filter 1789
    python3 collect_all_ground_truth.py --filter 1789 --redo
    python3 collect_all_ground_truth.py --input-dir /some/other/folder
"""

import sys
import csv
import argparse
from pathlib import Path

import numpy as np
import cv2

import matplotlib


def _select_backend():
    """
    Picks a GUI backend that is ACTUALLY USABLE on this machine.

    Note the subtlety that broke an earlier version: matplotlib.use()
    does not import the backend, it only records a preference. The
    real import happens later at plt.figure(), so wrapping use() in
    try/except catches nothing and the process dies at first figure
    with e.g. "ModuleNotFoundError: No module named 'tkinter'". So
    each candidate's underlying dependency is import-tested here,
    before committing to it.

    GTK3 is included as a fallback: it is present on this station and
    works fine once _patch_gtk3_keysym_bug() below is applied.
    """
    candidates = [
        ("TkAgg", "tkinter"),
        ("Qt5Agg", "PyQt5"),
        ("QtAgg", "PyQt5"),
        ("GTK3Agg", "gi"),
    ]
    import importlib
    for backend, dependency in candidates:
        try:
            importlib.import_module(dependency)
        except Exception:
            continue
        try:
            matplotlib.use(backend, force=True)
            return backend
        except Exception:
            continue
    return matplotlib.get_backend()


_ACTIVE_BACKEND = _select_backend()

import matplotlib.pyplot as plt          # noqa: E402  (must follow backend selection)
from matplotlib.widgets import Button    # noqa: E402


def _patch_gtk3_keysym_bug():
    """
    matplotlib's GTK3 backend calls
    cbook._unikey_or_keysym_to_mplkey(unikey, keysym) and that helper
    does keysym.lower() without a None check. Some key events (bare
    modifiers, certain Enter variants) arrive with keysym=None, so it
    raises inside matplotlib's own event handler -- which both spams
    the console and prevents the application's key callback from ever
    running. Wrap it so None is handled instead of raising.
    """
    try:
        from matplotlib import cbook
        original = cbook._unikey_or_keysym_to_mplkey

        def safe(unikey, keysym):
            try:
                return original(unikey, keysym)
            except AttributeError:
                return unikey if unikey else ""

        cbook._unikey_or_keysym_to_mplkey = safe
    except Exception:
        pass  # helper is private; if it moves, the buttons still work


_patch_gtk3_keysym_bug()


def verify_backend_usable():
    """
    Creates and immediately closes a throwaway figure to force the
    backend to actually import and connect to the display. Without
    this, an unusable backend only reveals itself as a long traceback
    at the first real image, after the user has already committed to
    a collection session.
    """
    try:
        fig = plt.figure(figsize=(1, 1))
        plt.close(fig)
        return True, None
    except Exception as e:
        return False, e


def report_backend_failure(error):
    print("=" * 78)
    print("ERROR: no usable interactive matplotlib backend on this machine.")
    print("=" * 78)
    print(f"Backend attempted : {_ACTIVE_BACKEND}")
    print(f"Underlying error  : {type(error).__name__}: {error}")
    print()
    print("Ground-truth collection is inherently point-and-click, so it needs a")
    print("working GUI backend and a display. Options, easiest first:")
    print()
    print("  1. Install Tk support (usually the quickest fix):")
    print("       sudo apt install python3-tk")
    print()
    print("  2. Install PyQt5 instead:")
    print("       pip3 install PyQt5")
    print()
    print("  3. If you are on SSH, make sure a display is available:")
    print("       ssh -X argus_user@<host>       # X forwarding")
    print("     or run this from the NUC's own keyboard/monitor.")
    print()
    print("  4. Copy the images to a machine that has a GUI, collect there,")
    print("     and copy the resulting CSVs back into ground_truth/.")
    print()
    print("Your images are already preserved in archive/images/, so nothing is")
    print("lost by deferring this -- collection can happen any time.")
    print("=" * 78)


DEFAULT_INPUT_FOLDER = "/mnt/I2Rgus_Data/waterline/input"
DEFAULT_GROUND_TRUTH_FOLDER = "/mnt/I2Rgus_Data/waterline/ground_truth"

# Mean grayscale brightness (0-255) below which an image is treated as
# a night frame with no visible waterline. Matches
# NIGHT_BRIGHTNESS_THRESHOLD in waterline_detector_v5.py. Night
# captures (e.g. 02:00 UTC = 22:00 local at Cape Cod) render as solid
# black, so there is nothing to click -- auto-skipping them avoids
# paging through blank windows one at a time.
NIGHT_BRIGHTNESS_THRESHOLD = 20.0


def measure_brightness(image):
    """Mean grayscale brightness (0-255) of a BGR image."""
    if image.ndim == 3:
        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    else:
        gray = image
    return float(np.mean(gray))


def collect_points_for_image(image_path, index, total,
                             display_scale=None, live_preview=True, dpi=70,
                             max_display_px=1100):
    """
    SPEED NOTES (this matters a lot over a slow/remote X connection):

      * The image is drawn ONCE. The previous version called
        ax.clear() + ax.imshow() on every click, which re-renders and
        re-transmits the entire full-resolution canvas across the X
        link per point -- the direct cause of multi-minute click
        latency over cellular. Now only a single line artist is
        updated.
      * The image is DOWNSCALED for display (max_display_px). Clicked
        coordinates are converted back to full-resolution before being
        saved, so accuracy of the stored ground truth is set by the
        original image, not the preview. At 2448px wide downscaled to
        1100px, that is ~5x fewer pixels to push per redraw.
      * live_preview=False skips canvas redraws on click entirely.
        Points are confirmed on the CONSOLE instead (a few bytes over
        SSH rather than megabytes over X11). You lose the on-image
        dots while clicking, but each click becomes effectively
        instant. The line is still interpolated and saved normally.

    Points are captured in display coordinates and converted to
    full-resolution coordinates on save.
    """
    image = cv2.imread(str(image_path))
    if image is None:
        print(f"  Could not load image: {image_path}, skipping.")
        return None, None, "skip"

    full_h, full_w = image.shape[:2]

    if display_scale is None:
        display_scale = min(1.0, float(max_display_px) / max(full_h, full_w))
    if display_scale < 1.0:
        disp = cv2.resize(
            image,
            (max(1, int(round(full_w * display_scale))),
             max(1, int(round(full_h * display_scale)))),
            interpolation=cv2.INTER_AREA,
        )
    else:
        display_scale = 1.0
        disp = image

    disp_rgb = cv2.cvtColor(disp, cv2.COLOR_BGR2RGB)
    disp_h, disp_w = disp_rgb.shape[:2]

    fig = plt.figure(figsize=(disp_w / dpi, (disp_h / dpi) + 1.1), dpi=dpi)
    ax = fig.add_axes([0.04, 0.13, 0.92, 0.82])

    points = []          # display coords
    action = {"result": None}

    ax.imshow(disp_rgb)  # drawn once, never re-drawn per click
    ax.set_xticks([])
    ax.set_yticks([])
    preview_line, = ax.plot([], [], "r-o", markersize=4, linewidth=1.5)

    mode = "live preview" if live_preview else "FAST (console feedback only)"
    ax.set_title(
        f"[{index}/{total}]  {image_path.name}\n"
        f"Left-click along the TRUE water/sand line, LEFT TO RIGHT.  "
        f"Right-click undoes.  Display {int(display_scale*100)}%  |  {mode}"
    )

    def refresh_line():
        if not live_preview:
            return
        if points:
            preview_line.set_data([p[0] for p in points], [p[1] for p in points])
        else:
            preview_line.set_data([], [])
        fig.canvas.draw_idle()

    def on_click(event):
        if event.inaxes != ax:
            return
        if event.button == 1:
            points.append((event.xdata, event.ydata))
            # Console feedback is ~40 bytes; an X11 canvas redraw is megabytes.
            print(f"    point {len(points):3d}  "
                  f"full-res x={event.xdata / display_scale:7.1f} "
                  f"y={event.ydata / display_scale:7.1f}")
            refresh_line()
        elif event.button == 3:
            if points:
                points.pop()
                print(f"    undo -> {len(points)} point(s)")
                refresh_line()

    def on_key(event):
        if event.key in ("enter", "return"):
            action["result"] = "done"
        elif event.key == "s":
            action["result"] = "skip"
        elif event.key == "q":
            action["result"] = "quit"

    def on_close(event):
        if action["result"] is None:
            action["result"] = "skip"

    fig.canvas.mpl_connect("button_press_event", on_click)
    fig.canvas.mpl_connect("key_press_event", on_key)
    fig.canvas.mpl_connect("close_event", on_close)

    def make_button(left, label, color, result):
        axb = fig.add_axes([left, 0.02, 0.19, 0.075])
        btn = Button(axb, label, color=color, hovercolor="#dddddd")

        def handler(_event):
            if result == "undo":
                if points:
                    points.pop()
                    print(f"    undo -> {len(points)} point(s)")
                    refresh_line()
                return
            action["result"] = result

        btn.on_clicked(handler)
        return btn

    buttons = [
        make_button(0.05, "Save & Next", "#b6e3b6", "done"),
        make_button(0.28, "Undo point", "#e8e8e8", "undo"),
        make_button(0.51, "Skip image", "#f0dca8", "skip"),
        make_button(0.74, "Quit", "#f0b8b8", "quit"),
    ]

    plt.show(block=False)

    while action["result"] is None:
        plt.pause(0.1)

    try:
        plt.close(fig)
    except Exception:
        pass

    del buttons

    # Convert display coords back to full-resolution coords.
    full_points = [(x / display_scale, y / display_scale) for x, y in points]
    return full_points, full_w, action["result"]


def build_full_width_line(points, width):
    """
    Interpolates the clicked points across every image column. Points
    arrive in FULL-RESOLUTION coordinates (collect_points_for_image
    converts them back from display coordinates before returning), so
    `width` must be the full-resolution image width -- the saved
    ground truth is therefore at native resolution regardless of what
    preview scale was used for clicking.
    """
    points = sorted(points, key=lambda p: p[0])
    xs = np.array([p[0] for p in points], dtype=np.float32)
    ys = np.array([p[1] for p in points], dtype=np.float32)
    all_columns = np.arange(width, dtype=np.float32)
    return np.interp(all_columns, xs, ys)


def save_ground_truth_csv(output_path, line):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Column", "Row"])
        for x, y in enumerate(line):
            writer.writerow([x, int(round(y))])


def main():
    parser = argparse.ArgumentParser(
        description="Collect ground-truth shorelines by clicking, with resume and filtering."
    )
    parser.add_argument("--input-dir", default=DEFAULT_INPUT_FOLDER)
    parser.add_argument("--ground-truth-dir", default=DEFAULT_GROUND_TRUTH_FOLDER)
    parser.add_argument("--filter", default=None,
                        help="Only collect for images whose filename contains this substring "
                             "(e.g. --filter 1789 for one day's captures).")
    parser.add_argument("--fast", action="store_true",
                        help="FAST MODE for slow/remote connections. Disables on-image "
                             "redraws while clicking; each point is confirmed on the console "
                             "instead. Removes the per-click canvas transfer that makes "
                             "clicking take minutes over X11 on cellular. Points are still "
                             "recorded and saved identically.")
    parser.add_argument("--max-display-px", type=int, default=1100,
                        help="Downscale the preview so its largest dimension is at most this "
                             "(default 1100). Clicks are converted back to full-resolution "
                             "coordinates, so saved ground truth keeps full accuracy. Lower "
                             "= less data over the wire = faster.")
    parser.add_argument("--dpi", type=int, default=70,
                        help="Figure DPI (default 70). Lower means fewer pixels to transmit.")
    parser.add_argument("--min-brightness", type=float, default=NIGHT_BRIGHTNESS_THRESHOLD,
                        help=f"Auto-skip images whose mean brightness is below this "
                             f"(default {NIGHT_BRIGHTNESS_THRESHOLD}). Night captures are solid "
                             f"black and have no clickable waterline. Set 0 to disable.")
    parser.add_argument("--redo", action="store_true",
                        help="Also re-collect images that already have a ground-truth CSV "
                             "(default is to skip them, so sessions can be resumed).")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    gt_dir = Path(args.ground_truth_dir)

    images = sorted(input_dir.glob("*.jpg"))
    if args.filter:
        images = [p for p in images if args.filter in p.name]

    if not images:
        print(f"No images found in {input_dir}"
              + (f" matching --filter {args.filter}" if args.filter else ""))
        sys.exit(1)

    already = []
    todo = []
    for p in images:
        if not args.redo and (gt_dir / (p.stem + ".csv")).exists():
            already.append(p)
        else:
            todo.append(p)

    dark = []
    if args.min_brightness > 0:
        kept = []
        for p in todo:
            img = cv2.imread(str(p))
            if img is None:
                kept.append(p)
                continue
            if measure_brightness(img) < args.min_brightness:
                dark.append(p)
            else:
                kept.append(p)
        todo = kept

    print(f"Backend in use: {_ACTIVE_BACKEND}")

    ok, backend_error = verify_backend_usable()
    if not ok:
        report_backend_failure(backend_error)
        sys.exit(1)

    print(f"Images matched     : {len(images)}")
    if already:
        print(f"Already collected  : {len(already)} (skipped -- use --redo to re-collect)")
    if dark:
        print(f"Night/dark frames  : {len(dark)} auto-skipped (mean brightness < "
              f"{args.min_brightness:g}; nothing to click)")
        for p in dark:
            print(f"                     - {p.name}")
    print(f"To collect now     : {len(todo)}")
    if args.fast:
        print("FAST MODE: on-image preview disabled; points confirm on the console.")
    print(f"Preview capped at  : {args.max_display_px}px, {args.dpi} dpi")
    print()

    if not todo:
        print("Nothing left to collect. Done.")
        return

    collected = 0
    skipped = 0

    for i, image_path in enumerate(todo, start=1):
        print(f"[{i}/{len(todo)}] {image_path.name}")

        points, width, result = collect_points_for_image(
            image_path, i, len(todo),
            live_preview=not args.fast,
            dpi=args.dpi,
            max_display_px=args.max_display_px,
        )

        if result == "quit":
            print("Quitting early -- everything collected so far is saved.")
            break

        if result == "skip" or points is None or len(points) < 2:
            if result != "skip":
                print("  Fewer than 2 points placed -- treated as skip.")
            else:
                print("  Skipped (no ground truth saved).")
            skipped += 1
            continue

        line = build_full_width_line(points, width)
        output_path = gt_dir / (image_path.stem + ".csv")
        save_ground_truth_csv(output_path, line)
        print(f"  Saved: {output_path.name} ({len(points)} points)")
        collected += 1

    print()
    print("=" * 60)
    print(f"Done. Collected {collected}, skipped {skipped}, out of {len(todo)} attempted.")
    if collected:
        print("Run archive_run.py now so this work survives the nightly wipe.")
    print("=" * 60)


if __name__ == "__main__":
    main()

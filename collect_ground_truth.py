#!/usr/bin/env python3
"""
Ground Truth Collection Tool (v2 - fixed finish/close handling)
-----------------------------------------------------------------
Click along the TRUE water/sand line in an image (left to right).

Controls:
    Left click   : add a point
    Right click  : undo last point
    ENTER key    : finish and save (press while the plot window is focused)
    Closing the window also finishes safely now.

Usage:
    python3 collect_ground_truth.py <image_path> <output_csv_path>
"""

import sys
import csv
import numpy as np
import cv2
import matplotlib.pyplot as plt
from pathlib import Path


def collect_points(image_path):
    image = cv2.imread(str(image_path))
    if image is None:
        raise RuntimeError(f"Could not load image: {image_path}")

    image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

    fig, ax = plt.subplots(figsize=(14, 9))
    ax.imshow(image_rgb)
    ax.set_title(
        "Click along the TRUE water/sand line, LEFT TO RIGHT.\n"
        "Right-click to undo. Press ENTER (window focused) or close window when done."
    )

    points = []
    finished = {"done": False}

    def redraw():
        ax.clear()
        ax.imshow(image_rgb)
        ax.set_title(
            "Click along the TRUE water/sand line, LEFT TO RIGHT.\n"
            "Right-click to undo. Press ENTER (window focused) or close window when done."
        )
        if points:
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            ax.plot(xs, ys, "r-o", markersize=4)
        fig.canvas.draw_idle()

    def on_click(event):
        if event.inaxes != ax:
            return
        if event.button == 1:  # left click
            points.append((event.xdata, event.ydata))
            redraw()
        elif event.button == 3:  # right click
            if points:
                points.pop()
                redraw()

    def on_key(event):
        if event.key == "enter":
            finished["done"] = True

    def on_close(event):
        finished["done"] = True

    fig.canvas.mpl_connect("button_press_event", on_click)
    fig.canvas.mpl_connect("key_press_event", on_key)
    fig.canvas.mpl_connect("close_event", on_close)

    print("Click points along the true shoreline, left to right.")
    print("Right-click to undo the last point.")
    print("Press ENTER (with the plot window focused) OR close the window when finished.")

    plt.show(block=False)

    while not finished["done"]:
        plt.pause(0.1)

    try:
        plt.close(fig)
    except Exception:
        pass

    if len(points) < 2:
        raise RuntimeError("Need at least 2 points to build a ground truth line.")

    width = image.shape[1]
    height = image.shape[0]

    return points, width, height


def build_full_width_line(points, width):
    points = sorted(points, key=lambda p: p[0])

    xs = np.array([p[0] for p in points], dtype=np.float32)
    ys = np.array([p[1] for p in points], dtype=np.float32)

    all_columns = np.arange(width, dtype=np.float32)
    interpolated = np.interp(all_columns, xs, ys)

    return interpolated


def save_ground_truth_csv(output_path, line):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Column", "Row"])
        for x, y in enumerate(line):
            writer.writerow([x, int(round(y))])

    print(f"Saved ground truth to: {output_path}")


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 collect_ground_truth.py <image_path> <output_csv_path>")
        sys.exit(1)

    image_path = sys.argv[1]
    output_csv_path = sys.argv[2]

    points, width, height = collect_points(image_path)

    print(f"Collected {len(points)} points.")

    line = build_full_width_line(points, width)

    save_ground_truth_csv(output_csv_path, line)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Plot Elevation Model
-----------------------
Overlays every captured shoreline for one camera on top of a real
background photo, color-coded by tide elevation at capture time. This
is a sanity-check "photographic model" of the elevation contour data
collected so far -- since georectification isn't done yet, this stays
in PIXEL space (like everything else in the project so far), but it's
still the fastest way to SEE how much of the beach face/tidal range has
actually been sampled, and where the gaps are, before investing in the
DEM-gridding step.

Usage:
    python3 plot_elevation_model.py <background_image> <contour_points_csv> <camera c1|c2> <output_png>

Example:
    python3 plot_elevation_model.py \\
        /mnt/I2Rgus_Data/waterline/input/1785155400.Mon.Jul.27_12_30_00.GMT.2026.CACO05.c1.snap.jpg \\
        /mnt/I2Rgus_Data/waterline/contour_points.csv \\
        c1 \\
        elevation_model_c1.png
"""

import sys
import csv
from collections import defaultdict

import numpy as np
import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize


def load_contour_points(path, camera):
    """
    Groups contour points by source_file (one group per capture
    instant), each with a constant tide_elevation and an array of
    (column, row) pixel points sorted by column.
    """
    frames = defaultdict(lambda: {"columns": [], "rows": [], "elevation": None})
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["camera"] != camera:
                continue
            key = row["source_file"]
            frames[key]["columns"].append(float(row["pixel_column"]))
            # float(), not int(): pixel_row carries sub-pixel precision
            # since detector v6.1. int() raised ValueError on "700.42".
            frames[key]["rows"].append(float(row["pixel_row"]))
            # Column name depends on the extractor version and source:
            # "tide_elevation_navd88" (current), or "tide_elevation"
            # (older output, before the datum was stated explicitly).
            if "tide_elevation_navd88" in row:
                frames[key]["elevation"] = float(row["tide_elevation_navd88"])
            elif "tide_elevation" in row:
                frames[key]["elevation"] = float(row["tide_elevation"])
            else:
                raise KeyError(
                    "No elevation column found in contour_points.csv. Expected "
                    "'tide_elevation_navd88' or 'tide_elevation'; got: "
                    + ", ".join(row.keys()))
            frames[key]["capture_time"] = row["capture_time_utc"]

    for key, data in frames.items():
        order = np.argsort(data["columns"])
        data["columns"] = np.array(data["columns"])[order]
        data["rows"] = np.array(data["rows"])[order]

    return frames


def main():
    if len(sys.argv) != 5:
        print("Usage: python3 plot_elevation_model.py <background_image> <contour_points_csv> <camera c1|c2> <output_png>")
        sys.exit(1)

    image_path, contour_csv, camera, output_path = sys.argv[1:5]

    image = cv2.imread(image_path)
    if image is None:
        print(f"ERROR: could not load background image: {image_path}")
        sys.exit(1)
    image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

    frames = load_contour_points(contour_csv, camera)
    if not frames:
        print(f"No contour points found for camera '{camera}' in {contour_csv}")
        sys.exit(1)

    elevations = [data["elevation"] for data in frames.values()]
    norm = Normalize(vmin=min(elevations), vmax=max(elevations))
    colormap = matplotlib.colormaps["turbo"]

    fig, ax = plt.subplots(figsize=(14, 9))
    ax.imshow(image_rgb)

    # Draw lowest-elevation shorelines first, highest last, so higher
    # (more "on top of the beach") lines aren't hidden under lower ones
    # when they happen to overlap in pixel space.
    for key, data in sorted(frames.items(), key=lambda kv: kv[1]["elevation"]):
        color = colormap(norm(data["elevation"]))
        ax.plot(data["columns"], data["rows"], "-", color=color, linewidth=2, alpha=0.85)

    ax.set_title(f"Camera {camera.upper()}: {len(frames)} captured shorelines, "
                 f"tide elevation {min(elevations):.2f}m to {max(elevations):.2f}m")
    ax.axis("off")

    sm = cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Tide elevation (m)")

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Plotted {len(frames)} shoreline(s) for camera '{camera}', "
          f"elevation range {min(elevations):.3f}m to {max(elevations):.3f}m")
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()

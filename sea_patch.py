#!/usr/bin/env python3
"""
A Fixed Patch Of Sea Surface, Seen By Any Camera Pointing
===========================================================
The optical wave gauge reads wave texture in the surf zone. To train on
2024-25 frames and run on today's camera, both must show the SAME patch of
sea at the SAME scale. Redrawing old frames into today's camera view only
does that if the two views overlap: with the lidar-corrected 2025
pointing, the Jan 2025 c2 camera saw just ~32% of today's view (47% of the
gauge's crop), so the earlier redraws put sea from the wrong place in most
of the frame.

A map view does not depend on the camera: every frame, old or live, is
projected onto one rectangle of sea surface -- a grid in metres aligned
with the beach (alongshore x cross-shore), at the frame's own water level.
Any pointing that saw the rectangle produces the same picture of it.

The rectangle is chosen once (choose_patch) as the largest one that every
pointing involved saw at every water level, and saved as JSON; training
(rectify_owg_images.py) and live inference (owg_live.py) both read it.

Image layout: row 0 is the seaward edge, so the sea is at the top as in the
camera; columns run alongshore in the direction of the shore bearing given
in the JSON (351 deg at Marconi: roughly north, to the right).
"""

import json
from pathlib import Path

import numpy as np

# Marconi, from the Jan 2025 lidar: the 0 m contour near N 4638400 and the
# beach trend (contours run ~351 deg, so seaward is ~81 deg).
DEFAULT_ORIGIN = (420150.0, 4638400.0)
DEFAULT_SHORE_BEARING = 351.0


class Patch:
    def __init__(self, origin, shore_bearing, along, cross, res):
        self.origin = np.array(origin, float)
        self.bearing = float(shore_bearing)
        self.along = (float(along[0]), float(along[1]))     # m, alongshore (u)
        self.cross = (float(cross[0]), float(cross[1]))     # m, seaward (v)
        self.res = float(res)
        b = np.radians(self.bearing)
        self.ua = np.array([np.sin(b), np.cos(b)])            # alongshore unit (E, N)
        self.va = np.array([np.sin(b + np.pi / 2), np.cos(b + np.pi / 2)])  # seaward unit

    @property
    def shape(self):
        return (int(round((self.cross[1] - self.cross[0]) / self.res)),
                int(round((self.along[1] - self.along[0]) / self.res)))

    def grid(self):
        """Ground E, N of every output pixel centre; row 0 = most seaward."""
        h, w = self.shape
        u = self.along[0] + (np.arange(w) + 0.5) * self.res
        v = self.cross[1] - (np.arange(h) + 0.5) * self.res
        U, V = np.meshgrid(u, v)
        E = self.origin[0] + U * self.ua[0] + V * self.va[0]
        N = self.origin[1] + U * self.ua[1] + V * self.va[1]
        return E, N

    def to_json(self):
        return {"origin": self.origin.round(3).tolist(), "shore_bearing_deg": self.bearing,
                "along_m": list(self.along), "cross_m": list(self.cross), "res_m": self.res,
                "shape_rows_cols": list(self.shape)}

    @classmethod
    def from_json(cls, d):
        if isinstance(d, (str, Path)):
            d = json.loads(Path(d).read_text())
        return cls(d["origin"], d["shore_bearing_deg"], d["along_m"], d["cross_m"], d["res_m"])

    def corners(self):
        out = []
        for u, v in ((self.along[0], self.cross[0]), (self.along[1], self.cross[0]),
                     (self.along[1], self.cross[1]), (self.along[0], self.cross[1])):
            out.append(self.origin + u * self.ua + v * self.va)
        return np.array(out)


def rectify(img, io, eo, z, patch):
    """Image -> patch at plane z. Returns (rectified image, fraction of the patch the camera saw)."""
    import cv2
    from view_reproject import ground_to_pixel
    E, N = patch.grid()
    U, V, ok = ground_to_pixel(E.ravel(), N.ravel(), z, io, eo)
    U = np.where(ok, U, -1).reshape(E.shape).astype(np.float32)
    V = np.where(ok, V, -1).reshape(E.shape).astype(np.float32)
    out = cv2.remap(img, U, V, cv2.INTER_LINEAR,
                    borderMode=cv2.BORDER_CONSTANT, borderValue=0)
    return out, float(ok.mean())


def largest_rectangle(mask):
    """Largest all-True axis-aligned rectangle: (r0, r1, c0, c1), end-exclusive."""
    h, w = mask.shape
    heights = np.zeros(w, dtype=int)
    best = (0, 0, 0, 0, 0)
    for r in range(h):
        heights = np.where(mask[r], heights + 1, 0)
        stack = []
        for c in range(w + 1):
            cur = heights[c] if c < w else 0
            start = c
            while stack and stack[-1][1] >= cur:
                s, hh = stack.pop()
                if hh * (c - s) > best[0]:
                    best = (hh * (c - s), r - hh + 1, r + 1, s, c)
                start = s
            stack.append((start, cur))
    return best[1:]


def largest_rectangle_aspect(mask, max_aspect=2.0):
    """
    Largest all-True axis-aligned rectangle whose sides differ by at most
    max_aspect: (r0, r1, c0, c1), end-exclusive. A thin strip has a large area
    but little wave texture in one direction, so area alone is the wrong goal.
    """
    h, w = mask.shape
    best = (0, 0, 0, 0, 0)
    for r0 in range(h):
        run_and = np.ones(w, bool)
        for r1 in range(r0 + 1, h + 1):
            run_and &= mask[r1 - 1]
            if not run_and.any():
                break
            rows = r1 - r0
            # longest run of True columns
            x = np.concatenate([[0], run_and.astype(np.int8), [0]])
            d = np.diff(x)
            starts, ends = np.where(d == 1)[0], np.where(d == -1)[0]
            k = np.argmax(ends - starts)
            width = min(ends[k] - starts[k], int(max_aspect * rows))
            if rows > max_aspect * width:
                continue
            if rows * width > best[0]:
                c0 = starts[k] + (ends[k] - starts[k] - width) // 2
                best = (rows * width, r0, r1, c0, c0 + width)
    return best[1:]


def choose_patch(views, origin=DEFAULT_ORIGIN, search=((-250.0, 250.0), (-60.0, 300.0)),
                 res=0.25, step=2.0, max_range=250.0, min_side=30.0, bearings=range(0, 180, 5),
                 max_aspect=2.0):
    """
    views: [(io, eo, [z, ...]), ...] -- every pointing and the water levels it is
    used at. Returns the largest rectangle (sides within max_aspect of each other,
    any orientation in `bearings`) that all of them see at every listed level,
    within max_range of each camera -- beyond that the pixel footprint is too
    coarse for wave texture -- at output resolution `res`.
    """
    from view_reproject import ground_to_pixel
    best = None
    for brg in bearings:
        probe = Patch(origin, brg, search[0], search[1], step)
        E, N = probe.grid()
        seen = np.ones(E.shape, bool)
        for io, eo, zs in views:
            rng = np.hypot(E - eo[0], N - eo[1])
            for z in zs:
                _, _, ok = ground_to_pixel(E.ravel(), N.ravel(), z, io, eo)
                seen &= ok.reshape(E.shape) & (rng <= max_range)
                if not seen.any():
                    break
        if not seen.any():
            continue
        r0, r1, c0, c1 = largest_rectangle_aspect(seen, max_aspect)
        area = (r1 - r0) * (c1 - c0) * step * step
        if best is None or area > best[0]:
            v_hi, v_lo = search[1][1] - r0 * step, search[1][1] - r1 * step
            u_lo, u_hi = search[0][0] + c0 * step, search[0][0] + c1 * step
            best = (area, Patch(origin, brg, (u_lo, u_hi), (v_lo, v_hi), res))
    if best is None:
        raise ValueError("no part of the search area is seen by every pointing")
    p = best[1]
    if min(p.along[1] - p.along[0], p.cross[1] - p.cross[0]) < min_side:
        raise ValueError(f"common area only {p.along[1] - p.along[0]:.0f} x "
                         f"{p.cross[1] - p.cross[0]:.0f} m")
    return p, best[0]

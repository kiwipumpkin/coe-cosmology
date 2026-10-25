#!/usr/bin/env python3
"""Match Chris's sketch: Time1 top / Time2 bottom, outer=S1, inner=S2, wavy tracks."""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "two_shells_time1_time2_vertical.png")


def wavy(ax, x0, y0, x1, y1, n=28, amp=0.055, color="0.35", lw=1.6):
    """Horizontal-ish wavy connector from (x0,y0) to (x1,y1) — but panels are vertical,
    so connect top→bottom with slight horizontal wiggle."""
    t = np.linspace(0, 1, n)
    xs = x0 + (x1 - x0) * t
    ys = y0 + (y1 - y0) * t
    # wiggle perpendicular to the main vertical direction
    dx = x1 - x0
    dy = y1 - y0
    L = np.hypot(dx, dy) + 1e-12
    px, py = -dy / L, dx / L
    w = amp * L * np.sin(2 * np.pi * 3.5 * t)
    xs = xs + w * px
    ys = ys + w * py
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="round", zorder=2, alpha=0.85)


def draw_shells(ax, cx, cy, R_in, R_out, D, scale_labels=1.0):
    """Draw CoE, S2(inner), S1(outer), US at cx+D. Chris labels: outer S1, inner S2."""
    # faint fill
    ax.add_patch(Circle((cx, cy), R_out, fill=True, fc="#4c72b0", alpha=0.05, zorder=0))
    ax.add_patch(Circle((cx, cy), R_in, fill=True, fc="white", zorder=1))
    ax.add_patch(Circle((cx, cy), R_in, fill=True, fc="#c44e52", alpha=0.05, zorder=1))

    ax.add_patch(Circle((cx, cy), R_in, fill=False, ec="#c44e52", lw=2.6, zorder=3))
    ax.add_patch(Circle((cx, cy), R_out, fill=False, ec="#4c72b0", lw=2.6, zorder=3))

    ax.plot(cx, cy, "o", color="0.1", ms=9, zorder=5)
    us_x = cx + D
    ax.plot(us_x, cy, "o", color="#2ca02c", ms=11, zorder=6)

    fs = 13 * scale_labels
    ax.text(cx - 0.35 * R_out, cy + 0.12 * R_out, "CoE", fontsize=fs, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.92), zorder=7)
    # S2 = inner (Chris), S1 = outer
    ax.text(cx, cy + R_in + 0.06 * R_out, "S2", fontsize=fs, ha="center", va="bottom", color="#c44e52",
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none", alpha=0.92), zorder=7)
    ax.text(cx, cy + R_out + 0.06 * R_out, "S1", fontsize=fs, ha="center", va="bottom", color="#4c72b0",
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none", alpha=0.92), zorder=7)
    ax.text(us_x, cy - 0.18 * R_out, "US", fontsize=fs, ha="center", va="top", color="#2ca02c",
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none", alpha=0.92), zorder=7)

    return {
        "coe": (cx, cy),
        "s2_top": (cx, cy + R_in),
        "s1_top": (cx, cy + R_out),
        "s2_right": (cx + R_in, cy),
        "s1_right": (cx + R_out, cy),
        "us": (us_x, cy),
        "s2_left": (cx - R_in, cy),
        "s1_left": (cx - R_out, cy),
    }


def main():
    # Same geometry as before; Chris: S1 outer, S2 inner
    R1_in, R1_out, D1 = 1.0, 2.4, 1.55
    s = 1.55
    R2_in, R2_out, D2 = R1_in * s, R1_out * s, D1 * s

    fig, ax = plt.subplots(figsize=(12, 17))
    ax.set_aspect("equal")
    ax.axis("off")

    # Place Time1 higher, Time2 lower — same CoE x so expansion is radial about shared vertical axis
    cx = 0.0
    cy1 = 5.8
    cy2 = -1.2

    p1 = draw_shells(ax, cx, cy1, R1_in, R1_out, D1)
    p2 = draw_shells(ax, cx, cy2, R2_in, R2_out, D2)

    # Wavy tracks Time1 → Time2 (like Chris sketch)
    wavy(ax, p1["s1_left"][0], p1["s1_left"][1], p2["s1_left"][0], p2["s1_left"][1], amp=0.04)
    wavy(ax, p1["s2_left"][0], p1["s2_left"][1], p2["s2_left"][0], p2["s2_left"][1], amp=0.05)
    wavy(ax, p1["us"][0], p1["us"][1] - 0.05, p2["us"][0], p2["us"][1] + 0.05, amp=0.045, color="#2ca02c")
    wavy(ax, p1["s2_right"][0], p1["s2_right"][1], p2["s2_right"][0], p2["s2_right"][1], amp=0.05)
    wavy(ax, p1["s1_right"][0], p1["s1_right"][1], p2["s1_right"][0], p2["s1_right"][1], amp=0.04)

    # Panel titles in open space
    ax.text(cx - R2_out - 0.3, cy1 + R1_out + 0.55, "Time 1", fontsize=16, ha="left", va="bottom",
            fontweight="bold")
    ax.text(cx - R2_out - 0.3, cy2 + R2_out + 0.55, "Time 2", fontsize=16, ha="left", va="bottom",
            fontweight="bold")

    ax.text(
        0.5, 0.01,
        "Same as your sketch, rotated: Time 1 above, Time 2 below · S1 outer, S2 inner · wavy lines = expansion tracks",
        transform=ax.transAxes, ha="center", va="bottom", fontsize=12, color="0.3",
    )
    ax.set_title(
        "Comoving shells: US stays between S2 and S1 as everything expands",
        fontsize=15, pad=18,
    )

    # margins
    ax.set_xlim(-R2_out - 1.2, R2_out + 1.4)
    ax.set_ylim(cy2 - R2_out - 0.9, cy1 + R1_out + 1.1)

    fig.tight_layout()
    fig.savefig(OUT, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", OUT)


if __name__ == "__main__":
    main()

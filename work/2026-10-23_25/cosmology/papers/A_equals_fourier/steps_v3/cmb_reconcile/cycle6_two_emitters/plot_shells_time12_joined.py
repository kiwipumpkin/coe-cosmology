#!/usr/bin/env python3
"""One continuous sketch: Time1 top → Time2 bottom joined by wavy expansion tracks (Chris style)."""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "two_shells_time1_time2_joined.png")


def wavy(ax, x0, y0, x1, y1, n=60, waves=4.0, amp=0.12, color="0.25", lw=1.8, z=4):
    t = np.linspace(0, 1, n)
    xs = x0 + (x1 - x0) * t
    ys = y0 + (y1 - y0) * t
    dx, dy = x1 - x0, y1 - y0
    L = np.hypot(dx, dy) + 1e-12
    px, py = -dy / L, dx / L
    w = amp * np.sin(2 * np.pi * waves * t)
    ax.plot(xs + w * px, ys + w * py, color=color, lw=lw, solid_capstyle="round",
            zorder=z, alpha=0.9)


def shells(ax, cx, cy, Rin, Rout, D, label_size=13):
    ax.add_patch(Circle((cx, cy), Rout, fill=False, ec="#4c72b0", lw=2.4, zorder=3))
    ax.add_patch(Circle((cx, cy), Rin, fill=False, ec="#c44e52", lw=2.4, zorder=3))
    ax.plot(cx, cy, "o", color="0.1", ms=8, zorder=5)
    usx = cx + D
    ax.plot(usx, cy, "o", color="#2ca02c", ms=10, zorder=6)
    ax.text(cx - 0.28 * Rout, cy + 0.1 * Rout, "CoE", fontsize=label_size, ha="right",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.95), zorder=7)
    ax.text(cx, cy + Rin + 0.05 * Rout, "S2", color="#c44e52", fontsize=label_size, ha="center",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.95), zorder=7)
    ax.text(cx, cy + Rout + 0.05 * Rout, "S1", color="#4c72b0", fontsize=label_size, ha="center",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.95), zorder=7)
    ax.text(usx + 0.08 * Rout, cy - 0.14 * Rout, "US", color="#2ca02c", fontsize=label_size,
            fontweight="bold", ha="left",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.95), zorder=7)
    return dict(
        coe=(cx, cy),
        us=(usx, cy),
        s2L=(cx - Rin, cy), s2R=(cx + Rin, cy),
        s1L=(cx - Rout, cy), s1R=(cx + Rout, cy),
        s2T=(cx, cy + Rin), s1T=(cx, cy + Rout),
    )


def main():
    # Time1 compact; Time2 expanded — share same CoE x, stack vertically with gap for tracks
    Rin1, Rout1, D1 = 1.0, 2.35, 1.55
    sc = 1.65
    Rin2, Rout2, D2 = Rin1 * sc, Rout1 * sc, D1 * sc

    cx = 0.0
    cy1 = 7.0
    # gap so wavy lines have room; Time2 CoE lower
    cy2 = cy1 - (Rout1 + Rout2) - 2.8

    fig, ax = plt.subplots(figsize=(11, 16))
    ax.set_aspect("equal")
    ax.axis("off")

    p1 = shells(ax, cx, cy1, Rin1, Rout1, D1)
    p2 = shells(ax, cx, cy2, Rin2, Rout2, D2)

    # Join corresponding points with wavy tracks (the whole point of the sketch)
    amp = 0.14
    wavy(ax, *p1["s1L"], *p2["s1L"], amp=amp, waves=3.5)
    wavy(ax, *p1["s2L"], *p2["s2L"], amp=amp * 0.9, waves=4.0)
    wavy(ax, *p1["us"], *p2["us"], amp=amp * 0.85, waves=3.8, color="#2ca02c", lw=2.0)
    wavy(ax, *p1["s2R"], *p2["s2R"], amp=amp * 0.9, waves=4.0)
    wavy(ax, *p1["s1R"], *p2["s1R"], amp=amp, waves=3.5)

    ax.text(cx - Rout2 - 0.6, cy1 + Rout1 + 0.35, "Time 1", fontsize=16, fontweight="bold")
    ax.text(cx - Rout2 - 0.6, cy2 + Rout2 + 0.35, "Time 2", fontsize=16, fontweight="bold")

    ax.set_title("One drawing: shells + US expand along the wavy tracks", fontsize=15, pad=12)
    ax.set_xlim(-Rout2 - 1.4, Rout2 + 1.6)
    ax.set_ylim(cy2 - Rout2 - 0.7, cy1 + Rout1 + 0.9)
    fig.tight_layout()
    fig.savefig(OUT, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close()
    print("wrote", OUT)

if __name__ == "__main__":
    main()

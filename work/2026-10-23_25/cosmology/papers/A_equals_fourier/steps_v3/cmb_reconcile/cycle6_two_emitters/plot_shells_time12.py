#!/usr/bin/env python3
"""Two-panel: concentric shells + US between at Time 1 and Time 2 (expanded, comoving)."""
from __future__ import annotations
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Arc

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "two_shells_time1_time2.png")


def draw_panel(ax, R_in, R_out, D, title, label_note):
    # CoE at origin
    ax.plot(0, 0, "o", color="0.15", ms=10, zorder=5)
    ax.text(-0.55 * R_out, 0.12 * R_out, "CoE", fontsize=14, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.92))

    # Shells
    c1 = Circle((0, 0), R_in, fill=False, ec="#c44e52", lw=2.8, zorder=3)
    c2 = Circle((0, 0), R_out, fill=False, ec="#4c72b0", lw=2.8, zorder=3)
    ax.add_patch(c1)
    ax.add_patch(c2)
    # soft fill between shells
    ring = Circle((0, 0), R_out, fill=True, fc="#4c72b0", alpha=0.06, zorder=0)
    hole = Circle((0, 0), R_in, fill=True, fc="white", alpha=1.0, zorder=1)
    ax.add_patch(ring)
    # redraw white disk then inner fill lightly
    inner_fill = Circle((0, 0), R_in, fill=True, fc="#c44e52", alpha=0.05, zorder=1)
    ax.add_patch(Circle((0, 0), R_in, fill=True, fc="white", zorder=1))
    ax.add_patch(inner_fill)

    # US between shells on +x
    ax.plot(D, 0, "o", color="#2ca02c", ms=12, zorder=6)
    ax.text(D, -0.22 * R_out, "US", fontsize=14, ha="center", va="top", color="#2ca02c",
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.92))

    # Labels for shells in open space
    ax.text(0, R_in + 0.08 * R_out, "S1", fontsize=14, ha="center", va="bottom", color="#c44e52",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9))
    ax.text(0, R_out + 0.08 * R_out, "S2", fontsize=14, ha="center", va="bottom", color="#4c72b0",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9))

    # Comoving outward arrows (same V feel — same length)
    arrow_len = 0.28 * R_out
    for r, yoff, col in [
        (R_in, 0.35 * R_in, "#c44e52"),
        (D, 0.0, "#2ca02c"),
        (R_out, -0.25 * R_out, "#4c72b0"),
    ]:
        # arrow along +x at height yoff near the object on the right side
        x0 = r * 0.92
        y0 = yoff if abs(yoff) > 1e-9 else 0.18 * R_out
        if abs(yoff) < 1e-9:
            y0 = 0.22 * R_out
        ax.annotate(
            "",
            xy=(x0 + arrow_len, y0),
            xytext=(x0, y0),
            arrowprops=dict(arrowstyle="-|>", color="0.45", lw=2.2, mutation_scale=14),
            zorder=4,
        )

    # Faint light rays toward US from both shells (along axis)
    # from S1 (left of US) and S2 (right of US)
    ax.annotate("", xy=(D - 0.08 * R_out, 0.05 * R_out), xytext=(R_in + 0.05 * R_out, 0.05 * R_out),
                arrowprops=dict(arrowstyle="->", color="#c44e52", lw=1.4, alpha=0.55), zorder=2)
    ax.annotate("", xy=(D + 0.08 * R_out, -0.05 * R_out), xytext=(R_out - 0.05 * R_out, -0.05 * R_out),
                arrowprops=dict(arrowstyle="->", color="#4c72b0", lw=1.4, alpha=0.55), zorder=2)

    ax.set_aspect("equal")
    pad = 0.35 * R_out
    ax.set_xlim(-R_out - pad, R_out + pad + 0.15 * R_out)
    ax.set_ylim(-R_out - pad, R_out + pad)
    ax.axis("off")
    ax.set_title(title, fontsize=16, pad=10, loc="left")
    ax.text(0.98, 0.02, label_note, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=12, color="0.35",
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="none", alpha=0.9))


def main():
    # Time 1 geometry (dimensionless)
    R1_in, R1_out, D1 = 1.0, 2.4, 1.55
    # Time 2: expanded ~1.55×, still US between (homologous-ish growth)
    s = 1.55
    R2_in, R2_out, D2 = R1_in * s, R1_out * s, D1 * s

    fig, axes = plt.subplots(2, 1, figsize=(14, 16))
    draw_panel(
        axes[0], R1_in, R1_out, D1,
        "Time 1 — concentric shells, US between them",
        "S1, US, S2 share ~same outward V",
    )
    draw_panel(
        axes[1], R2_in, R2_out, D2,
        "Time 2 — same layout, expanded together",
        "shells + US moved out; relative places unchanged",
    )
    fig.suptitle(
        "Comoving expansion: observer and both glowing shells coast outward together",
        fontsize=15, y=0.995,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(OUT, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", OUT)


if __name__ == "__main__":
    main()

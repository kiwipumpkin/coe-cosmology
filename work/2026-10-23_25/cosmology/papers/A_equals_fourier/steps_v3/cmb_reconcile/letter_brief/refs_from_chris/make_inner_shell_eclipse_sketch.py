#!/usr/bin/env python3
"""Inner-shell eclipse geometry sketch (cycle-6 concentric shells).

R_in=1, R_out=3, D=1.8 → α = arcsin(R_in/D) ≈ 33.7°.
Toward CoE, lines of sight with |θ|<α hit the inner shell only
(eclipsing the outer); outside that disk they hit the outer shell.

Regenerates cycle6_two_emitters/inner_shell_eclipse_sketch.png
to match the approved reference used in letter drafts.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc, Circle, FancyArrowPatch, Polygon, Wedge

HERE = Path(__file__).resolve().parent
OUT = HERE / "inner_shell_eclipse_sketch.png"

R_IN, R_OUT, D = 1.0, 3.0, 1.8
ALPHA = float(np.arcsin(R_IN / D))  # rad
ALPHA_DEG = float(np.rad2deg(ALPHA))

# Colors matching the approved sketch
C_SHELL = "#8B1A1A"
C_SHELL_FILL = "#F5E6E0"
C_US = "#1E5AA8"
C_ECLIPSE = "#B8C4D0"
C_CAPTION = "#333333"


def main() -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 8.2), dpi=150)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-3.65, 3.85)
    ax.set_ylim(-3.85, 4.55)

    # Title + params
    ax.text(
        0.0,
        4.28,
        "Inner shell eclipses outer toward CoE",
        ha="center",
        va="center",
        fontsize=14,
        fontweight="bold",
        color="0.05",
        clip_on=False,
    )
    ax.text(
        0.0,
        3.88,
        rf"$R_{{\rm in}}=1,\ R_{{\rm out}}=3,\ D=1.8$"
        rf"$\qquad\alpha=\arcsin(R_{{\rm in}}/D)\approx{ALPHA_DEG:.1f}^\circ$",
        ha="center",
        va="center",
        fontsize=11,
        color="0.15",
        clip_on=False,
    )

    # Outer fill + rings
    ax.add_patch(
        Circle((0, 0), R_OUT, facecolor=C_SHELL_FILL, edgecolor=C_SHELL, lw=2.2, zorder=1)
    )
    ax.add_patch(
        Circle((0, 0), R_IN, facecolor="white", edgecolor=C_SHELL, lw=2.2, zorder=2)
    )

    # Tangent geometry from US at (D, 0)
    # Tangent length and contact angle from CoE: sin(α)=R_in/D
    # Contact points on inner shell: angle ±(90°-α) from +x toward CoE side...
    # Observer at (D,0); tangents touch at angle φ where cos(φ_from_x_at_CoE related):
    # Contact: point = R_in * (cos ψ, sin ψ) with ψ = ±(π/2 + (π/2-α))? 
    # Standard: contact angles from +x (CoE→US): ψ = ±arccos(R_in/D) = ±(90°-α)... 
    # Actually: vector CoE→contact perpendicular to tangent; contact at angle ±β from US direction
    # where β = arcsin?  Contact angle from +x axis: ±(π - π/2 - α wait)
    # US at (D,0). Contact C on circle |C|=R_in, (C-US)·C = 0 ⇒ C·C = C·US ⇒ R_in^2 = D * C_x
    # ⇒ C_x = R_in^2/D, C_y = ±R_in*sqrt(1-(R_in/D)^2) = ±R_in*cos(α)
    cx = (R_IN**2) / D
    cy = R_IN * np.cos(ALPHA)
    c_top = np.array([cx, cy])
    c_bot = np.array([cx, -cy])
    us = np.array([D, 0.0])

    # Extend tangents past contact (toward CoE / beyond)
    def extend(p0, p1, scale=2.8):
        v = p1 - p0
        v = v / np.linalg.norm(v)
        return p0 + scale * v * R_OUT

    t_top_end = extend(us, c_top, scale=2.6)
    t_bot_end = extend(us, c_bot, scale=2.6)

    # Eclipse wedge (sector between tangents, toward CoE)
    # polygon: US → along top tangent far → arc or through contacts → bottom
    # Fill triangle-ish: US, far along top, far along bot (clipped visually by shells)
    # Redraw shells above eclipse fill for clean edges on near side
    ax.add_patch(
        Circle((0, 0), R_OUT, facecolor="none", edgecolor=C_SHELL, lw=2.2, zorder=3)
    )
    ax.add_patch(
        Circle((0, 0), R_IN, facecolor="white", edgecolor=C_SHELL, lw=2.2, zorder=4)
    )
    # re-fill eclipse only behind inner (clip): draw eclipse wedge again under labels
    # Keep simple: shade between tangents inside outer disk via wedge approx
    # Re-shade using Wedge from US is hard; use polygon clipped conceptually.

    # Tangents (dashed blue)
    for end in (t_top_end, t_bot_end):
        ax.plot(
            [us[0], end[0]],
            [us[1], end[1]],
            ls="--",
            color=C_US,
            lw=1.6,
            zorder=5,
        )

    # Direct LOS toward CoE (arrow to inner shell)
    ax.annotate(
        "",
        xy=(R_IN, 0),
        xytext=(D - 0.02, 0),
        arrowprops=dict(arrowstyle="-|>", color=C_US, lw=1.8, mutation_scale=12),
        zorder=6,
    )

    # Outer hit ray (solid red) — angle slightly outside α
    ang = ALPHA + np.deg2rad(18)
    # Ray from US in direction toward CoE-ish but above: direction = (-cos(ang), sin(ang))
    direction = np.array([-np.cos(ang), np.sin(ang)])
    # Intersect outer circle |US + t d|^2 = R_out^2
    # (us + t d)·(us + t d) = R_OUT^2
    a = np.dot(direction, direction)
    b = 2 * np.dot(us, direction)
    c = np.dot(us, us) - R_OUT**2
    disc = b**2 - 4 * a * c
    t1 = (-b - np.sqrt(disc)) / (2 * a)
    t2 = (-b + np.sqrt(disc)) / (2 * a)
    t_hit = t1 if t1 > 0 else t2
    hit = us + t_hit * direction
    ax.annotate(
        "",
        xy=hit,
        xytext=us,
        arrowprops=dict(arrowstyle="-|>", color=C_SHELL, lw=1.7, mutation_scale=11),
        zorder=6,
    )
    ax.plot(hit[0], hit[1], "o", color=C_SHELL, ms=5, zorder=7)
    ax.text(
        hit[0] + 0.28,
        hit[1] + 0.18,
        "hits outer",
        color=C_SHELL,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="bottom",
        zorder=9,
    )

    # Contact dots
    ax.plot(c_top[0], c_top[1], "o", color=C_SHELL, ms=5, zorder=7)
    ax.plot(c_bot[0], c_bot[1], "o", color=C_SHELL, ms=5, zorder=7)

    # Alpha arc
    arc_r = 0.55
    ax.add_patch(
        Arc(
            us,
            2 * arc_r,
            2 * arc_r,
            angle=0,
            theta1=180 - ALPHA_DEG,
            theta2=180,
            color=C_US,
            lw=1.5,
            zorder=8,
        )
    )
    ax.text(
        us[0] - arc_r * 0.85,
        0.22,
        r"$\alpha$",
        color=C_US,
        fontsize=13,
        ha="center",
        va="bottom",
        fontweight="bold",
    )

    # Labels
    ax.plot(0, 0, "o", color="k", ms=7, zorder=8)
    ax.text(0, -0.28, "CoE", ha="center", va="top", fontsize=12, fontweight="bold")
    ax.plot(D, 0, "o", color=C_US, ms=11, zorder=8)
    ax.text(D, -0.32, "us", color=C_US, ha="center", va="top", fontsize=13, fontweight="bold")

    ax.text(
        -1.55,
        R_OUT - 0.15,
        "outer shell",
        color=C_SHELL,
        ha="center",
        va="top",
        fontsize=11,
        fontweight="bold",
        zorder=9,
    )
    ax.text(
        0,
        R_IN + 0.12,
        "inner shell",
        color=C_SHELL,
        ha="center",
        va="bottom",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.text(
        -1.15,
        0.55,
        "hits inner\n(eclipses outer)",
        color=C_US,
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        zorder=9,
    )

    # Footer caption
    ax.text(
        0.0,
        -3.55,
        r"Sky disk $|\theta|<\alpha$ (and $360^\circ-\alpha$ twin): "
        r"lines of sight hit inner shell only",
        ha="center",
        va="center",
        fontsize=10.5,
        color=C_CAPTION,
        fontweight="bold",
    )

    fig.savefig(OUT, dpi=150, bbox_inches="tight", facecolor="white", pad_inches=0.15)
    plt.close(fig)
    print(f"wrote {OUT}  α={ALPHA_DEG:.2f}°")
    return OUT


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Two-emitter comoving figure matching letter att14_two_emitters look.

Title: "us between two emitters — comoving"
Top: parallel plates (comoving) with V (together)
Bottom: shells (comoving outward) with V_out ≈ same
Science axes / geometry unchanged in meaning (plates + shells, US between).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np

HERE = Path(__file__).resolve().parent
LETTER_FIGS = HERE.parent / "letter_brief" / "figs"
OUT_CYCLE = HERE / "two_emitters_comoving.png"
OUT_LETTER = LETTER_FIGS / "two_emitters_comoving.png"
OUT_DELIVERY = Path("/workspace/figs_letter_exact_python/Fig08_two_emitters.png")
GOLD = Path("/workspace/coe-figures/ref/att14_two_emitters.png")


def bulk_arrows(ax, y: float, x_positions: tuple[float, ...], length: float = 0.72) -> None:
    for x0 in x_positions:
        ax.annotate(
            "",
            xy=(x0 + length, y),
            xytext=(x0, y),
            arrowprops=dict(arrowstyle="-|>", color="0.40", lw=1.8, mutation_scale=11),
        )


def light_arrow(ax, x0: float, y: float, direction: int) -> None:
    span = 0.34
    if direction < 0:
        start, end = x0 - 0.06, x0 - span
    else:
        start, end = x0 + 0.06, x0 + span
    ax.annotate(
        "",
        xy=(end, y),
        xytext=(start, y),
        arrowprops=dict(arrowstyle="->", color="0.10", lw=1.05, mutation_scale=8),
    )


def draw_plates(ax) -> None:
    ax.set_xlim(-3.6, 3.6)
    ax.set_ylim(-2.45, 2.75)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("parallel plates (comoving)", fontsize=13, pad=8)

    bulk_arrows(ax, 2.20, (-2.45, -1.18, 0.09), length=0.72)
    ax.text(
        max((-2.45, -1.18, 0.09)) + 0.72 + 0.18,
        2.20,
        r"$V$ (together)",
        fontsize=11,
        color="0.30",
        ha="left",
        va="center",
    )

    y0, y1 = -1.55, 1.48
    x_left, x_right = -1.62, 1.62
    ax.plot([x_left, x_left], [y0, y1], color="0.05", lw=3.2, solid_capstyle="round")
    ax.plot([x_right, x_right], [y0, y1], color="0.05", lw=3.2, solid_capstyle="round")
    ax.text(x_left - 0.22, 0.0, "S2", fontsize=13, ha="right", va="center", fontweight="bold")
    ax.text(x_right + 0.22, 0.0, "S1", fontsize=13, ha="left", va="center", fontweight="bold")

    ax.plot(0.0, 0.0, "o", color="0.05", ms=10, zorder=5)
    ax.text(0.0, 0.26, "US", fontsize=12.5, ha="center", va="bottom", fontweight="bold")

    for x in (x_left, x_right):
        for y in np.linspace(-1.12, 1.05, 4):
            light_arrow(ax, x, y, -1)
            light_arrow(ax, x, y, +1)

    ax.text(
        0.5,
        -0.02,
        r"light both ways; plates + US share the same bulk $V$",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10.5,
        style="italic",
        color="0.25",
    )


def radial_light_arrow(ax, radius: float, angle_deg: float, direction: int) -> None:
    a = np.deg2rad(angle_deg)
    ux, uy = np.cos(a), np.sin(a)
    r0, r1 = radius + direction * 0.05, radius + direction * 0.34
    ax.annotate(
        "",
        xy=(r1 * ux, r1 * uy),
        xytext=(r0 * ux, r0 * uy),
        arrowprops=dict(arrowstyle="->", color="0.10", lw=1.0, mutation_scale=8),
    )


def draw_shells(ax) -> None:
    ax.set_xlim(-3.9, 4.4)
    ax.set_ylim(-3.25, 3.65)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("shells (comoving outward)", fontsize=13, pad=8)

    R1, R2 = 1.50, 2.65
    # Open on the left (west) so US mid-shell on the right stays visible
    theta1, theta2 = -110, 110
    for radius in (R1, R2):
        ax.add_patch(
            Arc(
                (0.0, 0.0),
                2 * radius,
                2 * radius,
                theta1=theta1,
                theta2=theta2,
                color="0.05",
                lw=2.9,
            )
        )

    ax.plot(0.0, 0.0, "+", color="0.10", ms=13, mew=1.9, zorder=5)
    ax.text(0.0, -0.36, "CoE", fontsize=12.5, ha="center", va="top")

    mid = 0.5 * (R1 + R2)
    ax.plot(mid, 0.0, "o", color="0.05", ms=10, zorder=6)
    ax.text(mid, 0.28, "US", fontsize=12.5, ha="center", va="bottom", fontweight="bold")

    for radius in (R1, R2):
        for angle in (-38, -17, 17, 38):
            radial_light_arrow(ax, radius, angle, -1)
            radial_light_arrow(ax, radius, angle, +1)

    a_lab = np.deg2rad(54)
    for radius, label in ((R1, "S1"), (R2, "S2")):
        ax.text(
            (radius + 0.14) * np.cos(a_lab),
            (radius + 0.14) * np.sin(a_lab),
            label,
            fontsize=12.5,
            ha="left",
            va="bottom",
            fontweight="bold",
        )

    for x0, y0, length in ((0.22, 0.58, 0.38), (0.92, 0.58, 0.38), (mid + 0.12, 0.0, 0.42)):
        ax.annotate(
            "",
            xy=(x0 + length, y0),
            xytext=(x0, y0),
            arrowprops=dict(arrowstyle="-|>", color="0.40", lw=1.7, mutation_scale=10),
        )
    ax.text(
        mid + 0.72,
        0.22,
        r"$V_{\mathrm{out}}\approx\mathrm{same}$",
        fontsize=11,
        color="0.30",
        ha="left",
        va="bottom",
    )

    ax.text(
        0.5,
        -0.02,
        r"shells + US coast outward together $\approx V_{\mathrm{rec}}$; relative Doppler $\approx 0$",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=10.5,
        style="italic",
        color="0.25",
    )


def main() -> None:
    # Target letter gold ~1181×1908 → figsize/dpi chosen to land near that after tight bbox
    fig = plt.figure(figsize=(7.2, 10.2), facecolor="white")
    fig.suptitle(
        "us between two emitters — comoving",
        fontsize=16,
        fontweight="bold",
        y=0.985,
    )
    gs = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=(1.0, 1.18),
        left=0.06,
        right=0.96,
        top=0.93,
        bottom=0.05,
        hspace=0.28,
    )

    ax_plates = fig.add_subplot(gs[0, 0])
    ax_shells = fig.add_subplot(gs[1, 0])
    draw_plates(ax_plates)
    draw_shells(ax_shells)

    outs = [OUT_CYCLE, OUT_LETTER, OUT_DELIVERY]
    for out in outs:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
        print("wrote", out)
    plt.close(fig)

    # Report dims vs gold
    from PIL import Image

    got = Image.open(OUT_DELIVERY).size
    gold = Image.open(GOLD).size if GOLD.exists() else None
    print("Fig08 size", got, "gold", gold)


if __name__ == "__main__":
    main()

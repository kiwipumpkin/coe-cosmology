#!/usr/bin/env python3
"""Donut-universe storyboard matching letter att14_storyboard (A–D).

Single stacked canvas (shared drawing scale) so fixed R_rec has the same
on-page size in every frame — matching att14. Vertical A→D.

Fill rule: S2 filled every frame (spherical plasma, inner red only).
Outer S1 = blue outline only — transparent, no disc / annulus fill.
Purple dashed = fixed recombination about CoE. Yellow dashed = visible at c about US.
Frame B must show S2 < R_recomb < S1 (only first shell past recombination).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

OUT = Path(__file__).resolve().parent / "donut_universe_storyboard.png"
OUT_DELIVERY = Path("/workspace/figs_letter_exact_python/Fig07_storyboard.png")
OUT_FIXED = Path("/workspace/figs_letter_exact_python/Fig07_storyboard_FIXED.png")
GOLD = Path("/workspace/coe-figures/ref/att14_storyboard.png")

C = 1.0
V_IN = 2.0
V_OUT = 4.0
V_D = 0.5 * (V_IN + V_OUT)
R_REC = 1.20
T_REC = R_REC / V_D
D_RECOMB = R_REC


def stage_at(t: float, title: str) -> dict:
    rin = V_IN * t
    rout = V_OUT * t
    d = V_D * t
    rvis = C * (t - T_REC) if d > D_RECOMB + 1e-12 else 0.0
    return dict(title=title, t=t, D=d, Rin=rin, Rout=rout, Rvis=rvis, Rrec=R_REC)


# Chris radii: A both inside R_rec; B with S2 < R_rec < S1; C/D both outside.
# Shared scale; fixed R_REC = 1.20. At t=0.40: Rin=0.80, Rout=1.60 → clear B ordering.
STAGES = [
    stage_at(0.12, "Time A — before recombination"),
    stage_at(0.40, "Time B — first shell past recombination"),
    stage_at(0.90, "Time C — both shells outside recombination"),
    stage_at(1.25, "Time D — today (donut universe)"),
]


def draw_stage_at(ax, s: dict, y0: float) -> float:
    """Draw one stage centred at (0, y0). Returns vertical half-extent used."""
    coe = (0.0, y0)
    us = (s["D"], y0)

    # S2 fill only
    ax.add_patch(
        Circle(coe, s["Rin"], facecolor="#f5c6c6", edgecolor="none", alpha=1.0, zorder=1)
    )
    ax.add_patch(Circle(coe, s["Rout"], fill=False, ec="#2980b9", lw=2.0, zorder=2))
    ax.add_patch(Circle(coe, s["Rin"], fill=False, ec="#c0392b", lw=1.8, zorder=2))
    ax.add_patch(
        Circle(coe, s["Rrec"], fill=False, ec="#8e44ad", lw=1.7, ls=(0, (5, 3)), zorder=3)
    )

    if s["Rvis"] > 1e-6:
        ax.add_patch(
            Circle(
                us,
                s["Rvis"],
                facecolor="#fef9e7",
                edgecolor="#f1c40f",
                lw=1.5,
                ls="--",
                alpha=0.40,
                zorder=3,
            )
        )
        ax.add_patch(Circle(us, s["Rvis"], fill=False, ec="#d4ac0d", lw=1.4, ls="--", zorder=4))

    ax.plot(*coe, "ko", ms=4.5, zorder=6)
    ax.plot(*us, "o", color="#27ae60", ms=6.0, zorder=6)

    # Labels
    ax.text(coe[0], coe[1] - 0.35, "CoE", fontsize=7.5, ha="center", va="top", color="#222")
    ax.text(
        us[0],
        us[1] + 0.32,
        "US",
        fontsize=8,
        ha="center",
        va="bottom",
        color="#1e8449",
        fontweight="bold",
    )

    if s["Rout"] > 0.35:
        ax.text(
            0.0,
            y0 + s["Rout"] * 0.85,
            "S1",
            fontsize=7.5,
            ha="center",
            va="bottom",
            color="#1a5276",
            fontweight="bold",
        )
        ax.text(
            -s["Rout"] * 0.72,
            y0 + s["Rout"] * 0.50,
            r"shells $>c$",
            fontsize=6.5,
            ha="center",
            va="center",
            color="0.45",
        )
    if s["Rin"] > 0.18:
        ax.text(
            0.0,
            y0 + max(s["Rin"] * 0.35, 0.12),
            "S2",
            fontsize=7.5,
            ha="center",
            va="center",
            color="#922b21",
            fontweight="bold",
        )

    ang = np.deg2rad(125)
    ax.text(
        s["Rrec"] * np.cos(ang),
        y0 + s["Rrec"] * np.sin(ang) + 0.05,
        "recomb. (fixed)",
        fontsize=7,
        color="#6c3483",
        ha="right",
        va="bottom",
    )

    if s["Rvis"] > 1e-6:
        ax.text(
            us[0] + max(s["Rvis"] * 0.30, 0.15),
            us[1] - s["Rvis"] * 0.80,
            "visible at $c$",
            fontsize=7,
            color="#b7950b",
            ha="left",
            va="top",
        )
        ax.text(
            us[0] + 0.12,
            us[1] + max(s["Rvis"] * 0.50, 0.18),
            "visible",
            fontsize=6.5,
            color="#b7950b",
            ha="left",
            va="bottom",
        )

    if s["Rin"] > 0.30:
        ax.annotate(
            "",
            xy=(-s["Rin"] - 0.40, y0 + 0.08),
            xytext=(-s["Rin"] + 0.02, y0 + 0.08),
            arrowprops=dict(arrowstyle="->", color="#7f8c8d", lw=1.3),
        )

    # Title just above this stage's bounding content
    top = max(s["Rout"], s["Rrec"], s["Rvis"] + 0.2, 0.9)
    ax.text(
        XMIN + 0.05,
        y0 + top + 0.35,
        s["title"],
        fontsize=10.5,
        fontweight="bold",
        ha="left",
        va="bottom",
        color="#111",
        clip_on=False,
    )
    return top


def main() -> None:
    # Lay out stages top→bottom with shared scale; gap between panels
    gap = 0.85
    title_pad = 0.55
    centers = []
    y_cursor = 0.0
    # First pass: compute vertical centers from top
    half_extents = []
    for s in STAGES:
        top = max(s["Rout"], s["Rrec"], s["Rvis"] + 0.15 if s["Rvis"] else 0.0, 0.85)
        half_extents.append(top)

    # Place from top: start with first title+content
    y = 0.0
    for i, (s, top) in enumerate(zip(STAGES, half_extents)):
        if i == 0:
            y = -(title_pad + top)
        else:
            prev_top = half_extents[i - 1]
            # previous center was centers[-1]; bottom of prev = centers[-1]-prev_top
            y = centers[-1] - prev_top - gap - title_pad - top
        centers.append(y)

    global XMIN, XMAX
    # Horizontal: fit largest stage + US + labels
    xmax = max(max(s["Rout"], s["D"] + max(s["Rvis"], 0.2)) for s in STAGES) + 0.70
    XMIN, XMAX = -xmax, xmax

    ymin = centers[-1] - half_extents[-1] - 0.40
    ymax = centers[0] + half_extents[0] + title_pad + 0.35

    # Gold 1045×1795
    fig = plt.figure(figsize=(5.225, 8.975), facecolor="white")
    ax = fig.add_axes([0.04, 0.038, 0.92, 0.915])
    fig.suptitle(
        r"Donut universe — recombination fixed; shells cross it ($>c$); US mid-shell; visible at $c$",
        fontsize=9.2,
        fontweight="bold",
        y=0.988,
    )

    for s, y0 in zip(STAGES, centers):
        draw_stage_at(ax, s, y0)

    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal")
    ax.axis("off")

    fig.text(
        0.5,
        0.008,
        r"S2 = plasma fill; S1 transparent. Purple = fixed recomb. Yellow = visible at $c$. US mid-shell. Shells $>c$.",
        ha="center",
        va="bottom",
        fontsize=7.8,
        color="0.25",
    )

    for out in (OUT, OUT_DELIVERY, OUT_FIXED):
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=200, facecolor="white")
        print("wrote", out)
    plt.close(fig)

    from PIL import Image

    got = Image.open(OUT_DELIVERY).size
    gold = Image.open(GOLD).size if GOLD.exists() else None
    print("Fig07 size", got, "gold", gold)
    print("ylim", ymin, ymax, "xlim", XMIN, XMAX)
    for s in STAGES:
        print(
            f"  {s['title'][:6]} t={s['t']:.2f} Rin={s['Rin']:.2f} Rout={s['Rout']:.2f} "
            f"D={s['D']:.2f} Rvis={s['Rvis']:.2f} Rrec={s['Rrec']:.2f} "
            f"both_in={s['Rout']<s['Rrec']} S2<Rrec<S1={s['Rin']<s['Rrec']<s['Rout']} "
            f"both_out={s['Rin']>s['Rrec']}"
        )


if __name__ == "__main__":
    main()

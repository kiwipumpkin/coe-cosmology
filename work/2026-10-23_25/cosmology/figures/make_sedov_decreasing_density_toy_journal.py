#!/usr/bin/env python3
"""
Journal-style two-panel figure: Chris's decreasing-density Sedov–Taylor toy universe.

Panel (a): geometry schematic (CoE, off-center observer, sky circle, labels).
Panel (b): exponential rho(r)/rho0 and normalized enclosed mass vs r/L.

Outputs PNG (300 dpi) + PDF vector + tiny JSON sidecar.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Circle

OUT = Path(__file__).resolve().parent
PNG = OUT / "sedov_decreasing_density_toy_journal.png"
PDF = OUT / "sedov_decreasing_density_toy_journal.pdf"
JSON_PATH = OUT / "sedov_decreasing_density_toy_journal.json"

# Print-friendly (Wong / Tol)
C_BLUE = "#0072B2"
C_ORANGE = "#E69F00"
C_GREEN = "#009E73"
C_VERM = "#D55E00"
C_SKY = "#56B4E9"
C_GRAY = "#555555"
C_BLACK = "#000000"

# Toy exponential: rho/rho0 = exp(-r/L); observer example D/L = 3
D_OVER_L = 3.0
DENSITY_FORMULA = "rho(r)/rho0 = exp(-r/L)"
MASS_FORMULA = "m(r)/M = 1 - exp(-X)*(1 + X + X^2/2), X = r/L"


def apply_journal_style() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 300,
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman", "CMU Serif", "DejaVu Serif", "Times"],
            "mathtext.fontset": "cm",
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.minor.visible": True,
            "ytick.minor.visible": True,
            "lines.linewidth": 1.6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def rho_over_rho0(X: np.ndarray) -> np.ndarray:
    return np.exp(-np.asarray(X, dtype=float))


def m_enclosed_fraction(X: np.ndarray) -> np.ndarray:
    """Normalized enclosed mass for rho ∝ exp(-r/L)."""
    X = np.asarray(X, dtype=float)
    return 1.0 - np.exp(-X) * (1.0 + X + 0.5 * X**2)


def panel_geometry(ax: plt.Axes) -> None:
    """(a) CoE--O--P geometry; same color/math theme as panel (b)."""
    ax.set_aspect("equal")
    ax.set_xlim(-0.35, 2.95)
    ax.set_ylim(-1.05, 1.45)
    ax.axis("off")

    D = 1.55
    R = 0.72
    obs = np.array([D, 0.0])
    coe = np.array([0.0, 0.0])

    ph_deg = 55.0
    ph = np.deg2rad(ph_deg)
    dhat = np.array([np.cos(ph), np.sin(ph)])
    src = obs + R * dhat
    rhat_P = src / np.linalg.norm(src)
    Dhat = np.array([1.0, 0.0])
    perp_d = np.array([-dhat[1], dhat[0]])
    perp_r = np.array([-rhat_P[1], rhat_P[0]])

    ax.add_patch(Circle(obs, R, fill=False, ec="#888888", lw=1.1, zorder=1))

    # D green (matches r=D in panel b); r dashed; d blue
    ax.plot([coe[0], obs[0]], [0.0, 0.0], color=C_GREEN, lw=1.5, zorder=2)
    ax.plot([coe[0], src[0]], [coe[1], src[1]], color="#666666", lw=1.1,
            ls="--", zorder=2)
    ax.plot([obs[0], src[0]], [obs[1], src[1]], color=C_BLUE, lw=1.6, zorder=3)

    ax.text(D * 0.50, -0.16, r"$D$", fontsize=10, color=C_GREEN, ha="center")
    mid_r = 0.45 * src
    ax.text(*(mid_r - 0.10 * perp_r), r"$r$", fontsize=10, color="#555555",
            ha="right", va="top")
    mid_d = obs + 0.55 * (src - obs)
    ax.text(*(mid_d + 0.08 * perp_d), r"$d$", fontsize=10, color=C_BLUE,
            ha="left", va="bottom")

    ax.plot(*coe, "o", color=C_VERM, ms=7, zorder=6,
            markeredgecolor="white", markeredgewidth=0.6)
    ax.text(-0.06, 0.14, r"CoE", fontsize=9.5, color=C_VERM,
            ha="right", va="bottom", fontweight="bold")
    ax.plot(*obs, marker="s", ms=7.5, mfc=C_BLUE, mec="white", mew=0.7,
            linestyle="None", zorder=7)
    ax.text(D - 0.02, 0.16, r"$O$", fontsize=11, color=C_BLUE,
            ha="right", va="bottom", fontweight="bold")
    ax.plot(*src, "o", color=C_ORANGE, ms=7, zorder=7,
            markeredgecolor="white", markeredgewidth=0.6)
    ax.text(src[0] + 0.10, src[1] - 0.12, r"$P$", fontsize=11,
            color=C_ORANGE, ha="left", va="top", fontweight="bold")

    ax.add_patch(Arc(obs, 0.40, 0.40, angle=0, theta1=0, theta2=ph_deg,
                     color=C_BLUE, lw=1.1, zorder=5))
    ax.text(D + 0.28, 0.08, r"$\varphi$", fontsize=11, color=C_BLUE,
            ha="center", va="bottom")

    def vector_arrow(start, end, color, lw=1.55):
        ax.annotate(
            "", xy=tuple(end), xytext=tuple(start),
            arrowprops=dict(
                arrowstyle="-|>", color=color, lw=lw,
                mutation_scale=12, shrinkA=0, shrinkB=0,
            ),
            zorder=8,
        )

    gap = 0.05
    vO_len, vP_len = 0.48, 0.46
    aO_len, aP_len = 0.30, 0.30

    vector_arrow(obs + gap * Dhat, obs + (gap + vO_len) * Dhat, C_BLUE)
    ax.text(
        *(obs + (gap + vO_len + 0.04) * Dhat + np.array([0.02, -0.16])),
        r"$\vec{v}_O = H\vec{D}$", fontsize=9, color=C_BLUE,
        ha="left", va="top",
    )
    vector_arrow(obs - gap * Dhat, obs - (gap + aO_len) * Dhat, C_VERM, lw=1.35)
    ax.text(
        *(obs - (gap + aO_len + 0.02) * Dhat + np.array([0.0, 0.14])),
        r"$\vec{a}_O$", fontsize=9, color=C_VERM, ha="right", va="bottom",
    )

    vector_arrow(src + gap * rhat_P, src + (gap + vP_len) * rhat_P, C_BLUE)
    ax.text(
        *(src + (gap + vP_len + 0.06) * rhat_P + 0.02 * perp_r),
        r"$\vec{v}_P = H\vec{r}$", fontsize=9, color=C_BLUE,
        ha="left", va="bottom",
    )
    vector_arrow(src - gap * rhat_P, src - (gap + aP_len) * rhat_P, C_VERM, lw=1.35)
    ax.text(
        *(src - (gap + aP_len + 0.08) * rhat_P - 0.06 * perp_r),
        r"$\vec{a}_P$", fontsize=9, color=C_VERM, ha="center", va="top",
    )

    ax.set_title("Geometry: off-center observer", fontsize=10, pad=6)
    ax.text(
        -0.02, 1.02, r"$\mathbf{(a)}$", transform=ax.transAxes,
        fontsize=12, fontweight="bold", va="bottom", ha="left",
    )



def panel_density_mass(ax: plt.Axes) -> None:
    """(b) Exponential density and enclosed mass; observer at D/L."""
    X = np.linspace(0.0, 8.0, 800)
    rho = rho_over_rho0(X)
    mfrac = m_enclosed_fraction(X)
    X_D = D_OVER_L

    ax2 = ax.twinx()

    (l1,) = ax.plot(
        X,
        rho,
        color=C_BLUE,
        lw=1.8,
        label=r"$\rho(r)/\rho_0 = e^{-r/L}$",
    )
    (l2,) = ax2.plot(
        X,
        mfrac,
        color=C_VERM,
        lw=1.8,
        label=r"$m(r)/M = 1-e^{-X}(1+X+X^{2}/2)$",
    )

    ax.axvline(X_D, color=C_GREEN, ls="--", lw=1.2, zorder=2)
    ax.plot(X_D, rho_over_rho0(X_D), "o", color=C_GREEN, ms=6, zorder=5)
    ax2.plot(X_D, m_enclosed_fraction(X_D), "s", color=C_GREEN, ms=6, zorder=5)

    ax.set_xlabel(r"$r / L$")
    ax.set_ylabel(r"$\rho(r)/\rho_0$", color=C_BLUE)
    ax2.set_ylabel(r"$m(r)/M$", color=C_VERM)
    ax.tick_params(axis="y", colors=C_BLUE, which="both")
    ax2.tick_params(axis="y", colors=C_VERM, which="both", direction="in")
    ax2.spines["right"].set_color(C_VERM)
    ax.spines["left"].set_color(C_BLUE)

    ax.set_xlim(0, 8)
    ax.set_ylim(0, 1.05)
    ax2.set_ylim(0, 1.05)
    ax.grid(False)

    # Thin black spines on bottom/top; colored left/right already set
    for side in ("top", "bottom"):
        ax.spines[side].set_color(C_BLACK)
        ax.spines[side].set_linewidth(0.8)

    lines = [l1, l2, Line2D([0], [0], color=C_GREEN, ls="--", lw=1.2)]
    labels = [
        l1.get_label(),
        l2.get_label(),
        rf"$r=D$ ($D/L={D_OVER_L:g}$)",
    ]
    ax.legend(
        lines,
        labels,
        loc="center right",
        frameon=True,
        fancybox=False,
        framealpha=0.95,
        edgecolor="#bbbbbb",
        fontsize=7.5,
        borderpad=0.4,
        handlelength=1.8,
        labelspacing=0.4,
    )

    ax.set_title("Density and enclosed mass", fontsize=10, pad=6)
    ax.text(
        -0.02,
        1.02,
        r"$\mathbf{(b)}$",
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def main() -> None:
    apply_journal_style()

    # ~7 inch wide, two panels side by side
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.35), constrained_layout=False)
    fig.patch.set_facecolor("white")

    panel_geometry(axes[0])
    panel_density_mass(axes[1])

    fig.subplots_adjust(left=0.02, right=0.94, bottom=0.14, top=0.88, wspace=0.32)

    fig.savefig(PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(PDF, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    meta = {
        "png": str(PNG),
        "pdf": str(PDF),
        "script": str(OUT / "make_sedov_decreasing_density_toy_journal.py"),
        "density_formula": DENSITY_FORMULA,
        "enclosed_mass_formula": MASS_FORMULA,
        "D_over_L": D_OVER_L,
        "note": "Toy exponential panel uses D/L=3 (not recommended_profile twoscale x=5.9).",
        "dpi_png": 300,
        "figsize_inch": [7.0, 3.35],
    }
    JSON_PATH.write_text(json.dumps(meta, indent=2) + "\n")
    print(f"Wrote {PNG}")
    print(f"Wrote {PDF}")
    print(f"Wrote {JSON_PATH}")


if __name__ == "__main__":
    main()

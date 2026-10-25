#!/usr/bin/env python3
"""
Illustrative figures for Chris's Sedov–Taylor centered-universe model.

Binding physics:
  - CoE beyond horizon; observer off-center at D
  - Decreasing ρ(r) (Phase A eq. 3.13 exponential); Gauss g = -G M(<r)/r²
  - Pure homologous ST → isotropic H_obs; gravity → anisotropy
  - At d ≪ D: leading quadrupole P2; κ_natural = 0.15
  - Goal: first-principles H_obs(n̂) / A(l) vs Migkas angular data — NO A cosφ dipole overlay

Outputs PNGs + captions.md under this directory. Optional combined PDF.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Arc
from matplotlib.backends.backend_pdf import PdfPages

OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)

# Colorblind-friendly palette (Wong / Tol-ish)
C_BLUE = "#0072B2"
C_ORANGE = "#E69F00"
C_GREEN = "#009E73"
C_VERM = "#D55E00"
C_PURPLE = "#CC79A7"
C_SKY = "#56B4E9"
C_BLACK = "#000000"
C_GRAY = "#666666"

KAPPA_NAT = 0.15
ALPHA = 1.0
EPS = 0.1  # d/D


def style():
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 200,
            "font.size": 13,
            "axes.labelsize": 14,
            "axes.titlesize": 15,
            "legend.fontsize": 11,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.linewidth": 1.2,
            "lines.linewidth": 2.2,
            "mathtext.fontset": "dejavusans",
        }
    )


def m_enclosed_fraction(X: np.ndarray) -> np.ndarray:
    """Eq 3.13: m(r)/M = 1 - e^{-X}(1 + X + X²/2), X = r/r_s."""
    X = np.asarray(X, dtype=float)
    return 1.0 - np.exp(-X) * (1.0 + X + 0.5 * X**2)


def rho_over_rho0(X: np.ndarray) -> np.ndarray:
    """ρ ∝ e^{-r/r_s} ⇒ ρ/ρ0 = e^{-X}."""
    return np.exp(-np.asarray(X, dtype=float))


def frac_tidal(mu: np.ndarray, kappa: float, alpha: float, eps: float) -> np.ndarray:
    """
    Series for H_obs/H − 1 (point-mass / exterior tidal limit).
    Leading: ακ(3μ²−1) = 2ακ P2(μ)
    O(ε): ακ ε (4.5 μ − 7.5 μ³) → projects mainly on P3
    """
    mu = np.asarray(mu, dtype=float)
    lead = alpha * kappa * (3.0 * mu**2 - 1.0)
    odd = alpha * kappa * eps * (4.5 * mu - 7.5 * mu**3)
    return lead + odd


def half_ptp(kappa: float, alpha: float = ALPHA, eps: float = EPS, n: int = 4001) -> float:
    mu = np.linspace(-1.0, 1.0, n)
    f = frac_tidal(mu, kappa, alpha, eps)
    return 0.5 * (f.max() - f.min())


# ---------------------------------------------------------------------------
# Fig 1 — geometry
# ---------------------------------------------------------------------------
def fig_geometry(path: Path) -> None:
    """Improved toy-universe schematic (publication-style)."""
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(9.2, 7.0))
    ax.set_aspect("equal")
    ax.set_xlim(-0.55, 3.05)
    ax.set_ylim(-1.55, 1.75)
    ax.axis("off")

    D = 1.70
    R = 0.78
    obs = np.array([D, 0.0])
    coe = np.array([0.0, 0.0])

    radii = np.linspace(2.6, 0.12, 28)
    for rad in radii:
        t = 1.0 - rad / radii[0]
        a = 0.02 + 0.22 * (t ** 1.4)
        ax.add_patch(Circle(coe, rad, facecolor=C_SKY, edgecolor="none", alpha=a, zorder=0))

    R_hor = 1.35
    ax.add_patch(Circle(obs, R_hor, fill=False, ec=C_GRAY, lw=1.3, ls=(0, (4, 3)), zorder=1, alpha=0.85))
    ax.text(D + 0.15, R_hor - 0.08, "schematic horizon", fontsize=9.5, color=C_GRAY, style="italic", ha="left")

    for ang in np.linspace(0, 2 * np.pi, 12, endpoint=False):
        r0, r1 = 0.55, 0.78
        ax.annotate(
            "",
            xy=(r1 * np.cos(ang), r1 * np.sin(ang)),
            xytext=(r0 * np.cos(ang), r0 * np.sin(ang)),
            arrowprops=dict(arrowstyle="-|>", color=C_SKY, lw=1.0, mutation_scale=8, alpha=0.55),
            zorder=1,
        )

    ax.plot(*coe, "o", color=C_VERM, ms=14, zorder=6, markeredgecolor="white", markeredgewidth=1.2)
    ax.text(-0.12, 0.22, "CoE", fontsize=15, fontweight="bold", color=C_VERM, ha="right")
    ax.text(-0.12, 0.05, "(beyond horizon)", fontsize=9, color=C_VERM, ha="right", style="italic")

    ax.plot(*obs, "s", color=C_BLUE, ms=13, zorder=6, markeredgecolor="white", markeredgewidth=1.0)
    ax.text(D + 0.08, 0.18, "Observer", fontsize=14, fontweight="bold", color=C_BLUE, ha="left")
    ax.text(D + 0.08, 0.02, "(us)", fontsize=10, color=C_BLUE, ha="left", style="italic")

    ax.annotate(
        "",
        xy=(D - 0.12, 0),
        xytext=(0.18, 0),
        arrowprops=dict(arrowstyle="-|>", color=C_GREEN, lw=2.6, mutation_scale=16),
        zorder=3,
    )
    ax.text(D * 0.48, -0.22, r"$\hat{D}$", fontsize=16, color=C_GREEN, ha="center", fontweight="bold")
    ax.text(D * 0.48, -0.40, r"$D$", fontsize=15, color=C_GREEN, ha="center")

    ax.annotate(
        "",
        xy=(D - 0.62, 0.0),
        xytext=(D - 0.18, 0.0),
        arrowprops=dict(arrowstyle="-|>", color=C_PURPLE, lw=1.6, mutation_scale=11),
        zorder=4,
    )
    ax.text(D - 0.95, 0.32, r"toward CoE", fontsize=10.5, color=C_PURPLE)
    ax.text(D - 0.95, 0.16, r"$(\varphi=\pi)$", fontsize=9.5, color=C_PURPLE)

    ax.annotate(
        "",
        xy=(D + 0.55, 0.0),
        xytext=(D + 0.12, 0.0),
        arrowprops=dict(arrowstyle="-|>", color=C_GRAY, lw=1.3, mutation_scale=10),
        zorder=4,
    )
    ax.text(D + 0.58, -0.28, r"anti-center $(\varphi=0)$", fontsize=9.5, color=C_GRAY, ha="left")

    ax.add_patch(Circle(obs, R, fill=False, ec=C_BLACK, lw=2.0, zorder=2))

    phis = np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315])
    for ph in phis:
        s = obs + R * np.array([np.cos(ph), np.sin(ph)])
        ax.plot(*s, "o", color=C_ORANGE, ms=8, zorder=5, markeredgecolor="white", markeredgewidth=0.6)

    ph_deg = 55.0
    ph = np.deg2rad(ph_deg)
    src = obs + R * np.array([np.cos(ph), np.sin(ph)])
    ax.plot(*src, "o", color=C_VERM, ms=11, zorder=7, markeredgecolor="white", markeredgewidth=1.0)
    ax.text(src[0] + 0.06, src[1] + 0.08, "source", fontsize=10, color=C_VERM, fontweight="bold")

    ax.plot([obs[0], src[0]], [obs[1], src[1]], color=C_ORANGE, lw=2.2, zorder=4)
    ax.plot([coe[0], src[0]], [coe[1], src[1]], color=C_GRAY, lw=1.5, ls="--", zorder=3)

    mid_d = obs + 0.55 * (src - obs)
    ax.text(mid_d[0] + 0.10, mid_d[1] + 0.02, r"$d$", fontsize=15, color=C_ORANGE, fontweight="bold")
    mid_r = 0.42 * src
    ax.text(mid_r[0] - 0.18, mid_r[1] + 0.14, r"$r$", fontsize=15, color=C_GRAY, fontweight="bold")

    ax.add_patch(Arc(obs, 0.62, 0.62, angle=0, theta1=0, theta2=ph_deg, color=C_BLACK, lw=1.6, zorder=5))
    ax.text(D + 0.42, 0.16, r"$\varphi$", fontsize=16, fontweight="bold")

    ax.text(D + 0.35, -R * 0.92, r"$R\sim d$", fontsize=12, color=C_BLACK)

    ax.annotate(
        r"decreasing $\rho(r)$" + "\n" + r"(centrally concentrated)",
        xy=(0.55, -0.95),
        xytext=(0.85, -1.35),
        fontsize=10.5,
        color=C_SKY,
        ha="left",
        arrowprops=dict(arrowstyle="-|>", color=C_SKY, lw=1.2, mutation_scale=10),
    )

    legend_elems = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_VERM, markersize=10, label="CoE (origin)"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor=C_BLUE, markersize=9, label="Observer at $(D,0)$"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=C_ORANGE, markersize=8, label="Sources on sky circle"),
        Line2D([0], [0], color=C_GREEN, lw=2.5, label=r"Preferred axis $\hat{D}$"),
        Line2D([0], [0], color=C_GRAY, lw=1.5, ls="--", label=r"CoE–source $r$"),
        Line2D([0], [0], color=C_ORANGE, lw=2.0, label=r"Observer–source $d$"),
    ]
    ax.legend(
        handles=legend_elems,
        loc="upper right",
        frameon=True,
        fancybox=True,
        framealpha=0.95,
        edgecolor="#cccccc",
        fontsize=9.5,
        borderpad=0.6,
    ).set_zorder(10)

    ax.set_title(
        "Toy centered Universe: off-center observer in a Sedov–Taylor blast about the CoE",
        fontsize=13,
        pad=10,
        color=C_BLACK,
    )
    fig.text(
        0.5,
        0.02,
        r"Pure homology $\vec{v}=H\vec{r}$ $\Rightarrow$ isotropic $H_{\mathrm{obs}}$; Gauss gravity with decreasing $\rho$ breaks homology.",
        ha="center",
        fontsize=9.5,
        color=C_GRAY,
        style="italic",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


# Fig 2 — density & mass
# ---------------------------------------------------------------------------
def fig_density_mass(path: Path) -> None:
    X = np.linspace(0, 8, 600)
    rho = rho_over_rho0(X)
    mfrac = m_enclosed_fraction(X)
    # Example: mark r = D. Choose D/r_s ≈ 3 as a typical “soft exterior” example.
    X_D = 3.0

    fig, ax1 = plt.subplots(figsize=(8.0, 5.2))
    ax2 = ax1.twinx()

    (l1,) = ax1.plot(X, rho, color=C_BLUE, label=r"$\rho(r)/\rho_0 = e^{-r/r_s}$")
    (l2,) = ax2.plot(
        X,
        mfrac,
        color=C_VERM,
        label=r"$m(r)/M = 1-e^{-X}(1+X+X^2/2)$",
    )

    ax1.axvline(X_D, color=C_GREEN, ls="--", lw=1.8, label=r"$r=D$ (example, $D/r_s=3$)")
    ax1.plot(X_D, rho_over_rho0(X_D), "o", color=C_GREEN, ms=8, zorder=5)
    ax2.plot(X_D, m_enclosed_fraction(X_D), "s", color=C_GREEN, ms=8, zorder=5)

    ax1.set_xlabel(r"$r / r_s$")
    ax1.set_ylabel(r"$\rho(r)/\rho_0$", color=C_BLUE)
    ax2.set_ylabel(r"$m(r)/M$", color=C_VERM)
    ax1.tick_params(axis="y", colors=C_BLUE)
    ax2.tick_params(axis="y", colors=C_VERM)
    ax1.set_xlim(0, 8)
    ax1.set_ylim(0, 1.05)
    ax2.set_ylim(0, 1.05)
    ax1.grid(True, alpha=0.3)

    lines = [l1, l2]
    labels = [ln.get_label() for ln in lines]
    # add vline to legend
    from matplotlib.lines import Line2D

    lines.append(Line2D([0], [0], color=C_GREEN, ls="--", lw=1.8))
    labels.append(r"$r=D$ example ($D/r_s=3$)")
    ax1.legend(lines, labels, loc="center right", framealpha=0.95)

    ax1.set_title("Phase A eq. (3.13): exponential decreasing density & enclosed mass")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 3 — H_obs vs angle
# ---------------------------------------------------------------------------
def load_migkas_fig6():
    """Digitised Migkas+2020 Fig.6 A/Aall and B/Ball vs Galactic longitude."""
    dig = Path("/workspace/cosmology/migkas/digitised")
    A = np.genfromtxt(
        dig / "migkas2020_fig6_A_vs_longitude.csv",
        delimiter=",",
        names=True,
    )
    B = np.genfromtxt(
        dig / "migkas2020_fig6_B_vs_longitude.csv",
        delimiter=",",
        names=True,
    )
    return A, B


def model_A_over_Aall_vs_longitude(
    l_deg: np.ndarray,
    l_axis_deg: float = 280.0,
    kappa: float = KAPPA_NAT,
    alpha: float = ALPHA,
    eps: float = EPS,
) -> np.ndarray:
    """
    Map first-principles H_obs anisotropy onto Migkas A/Aall vs Galactic longitude.

    Preferred axis on Galactic plane at longitude l_axis (default 280 deg =
    Migkas low-A / high-H0 direction). At b=0: mu = cos(l - l_axis).

    LX-T: A * H0^2 conserved => A/A_iso = 1/(1+f)^2 with f = H_obs/H - 1.
    Renormalise so mean_l(A/Aall) = 1 to match Fig.6 definition.
    """
    l = np.asarray(l_deg, dtype=float)
    mu = np.cos(np.deg2rad(l - l_axis_deg))
    f = frac_tidal(mu, kappa, alpha, eps)
    A_raw = 1.0 / (1.0 + f) ** 2
    return A_raw / np.mean(A_raw)


def fig_Hobs_vs_angle(path: Path) -> None:
    """Migkas Fig.6 digitised A,B + model on Galactic longitude. No dipole overlay."""
    Adata, Bdata = load_migkas_fig6()
    l_grid = np.linspace(0.0, 360.0, 721)
    A_model = model_A_over_Aall_vs_longitude(l_grid, l_axis_deg=280.0)
    # Fair comparison: Migkas Fig.6 uses Delta-l = 90 deg sliding bins
    dl = l_grid[1] - l_grid[0]
    half = int(round(45.0 / dl))
    kern = np.ones(2 * half + 1)
    kern /= kern.sum()
    A_smooth = np.convolve(np.pad(A_model, (half, half), mode="wrap"), kern, mode="valid")

    fig, (axA, axB) = plt.subplots(
        2, 1, figsize=(8.4, 7.2), sharex=True, gridspec_kw={"height_ratios": [1.35, 1.0]}
    )

    axA.axhline(
        1.0,
        color=C_GREEN,
        lw=1.6,
        label=r"(a) $M=0$ homology $\Rightarrow A/A_{\mathrm{all}}=1$",
    )
    axA.plot(
        l_grid,
        A_model,
        color=C_SKY,
        lw=1.5,
        alpha=0.75,
        label=r"(b) model unsmoothed ($b=0$ cut)",
    )
    axA.plot(
        l_grid,
        A_smooth,
        color=C_BLUE,
        lw=2.6,
        label=rf"(c) model $\times\Delta\ell=90^\circ$ bin (Fig.6-like)"
        + "\n"
        + rf"($\kappa={KAPPA_NAT}$, axis $\ell=280^\circ$)",
    )
    axA.errorbar(
        Adata["galactic_longitude_deg"],
        Adata["A_over_Aall"],
        yerr=Adata["A_over_Aall_err_1sigma"],
        fmt="o",
        color=C_VERM,
        ms=5,
        capsize=2,
        lw=1.0,
        zorder=5,
        label=r"Migkas+2020 Fig.6 $A/A_{\mathrm{all}}$ (digitised)",
    )
    axA.axvline(280, color=C_GRAY, ls=":", lw=1.2, alpha=0.8)
    axA.text(
        282,
        1.22,
        "low-A / high-H0\naxis ~280 deg",
        fontsize=10,
        color=C_GRAY,
    )
    axA.set_ylabel(r"$A/A_{\mathrm{all}}$")
    axA.set_ylim(0.70, 1.30)
    axA.legend(loc="upper right", framealpha=0.95, fontsize=10)
    axA.grid(True, alpha=0.3)
    axA.set_title(
        r"Angular match: first-principles model vs Migkas $A(\ell)$ — no dipole template"
    )

    axB.axhline(
        1.0, color=C_GREEN, lw=1.6, label=r"$B/B_{\mathrm{all}}=1$ (isotropic slope)"
    )
    axB.errorbar(
        Bdata["galactic_longitude_deg"],
        Bdata["B_over_Ball"],
        yerr=Bdata["B_over_Ball_err_1sigma"],
        fmt="s",
        color=C_ORANGE,
        ms=5,
        capsize=2,
        lw=1.0,
        zorder=5,
        label=r"Migkas+2020 Fig.6 $B/B_{\mathrm{all}}$ (digitised)",
    )
    axB.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    axB.set_ylabel(r"$B/B_{\mathrm{all}}$")
    axB.set_xlim(0, 360)
    axB.set_xticks([0, 90, 180, 270, 360])
    axB.set_ylim(0.70, 1.30)
    axB.legend(loc="upper right", framealpha=0.95, fontsize=10)
    axB.grid(True, alpha=0.3)
    axB.text(
        0.02,
        0.05,
        r"Model does not predict $B(\ell)$; panel shows Migkas slope isotropy check only.",
        transform=axB.transAxes,
        fontsize=10,
        style="italic",
        color=C_GRAY,
    )

    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_model_vs_migkas_schematic(path: Path) -> None:
    """H_obs/H-1 vs angle from axis; Migkas delta-H from digitised A. Qualitative."""
    Adata, _ = load_migkas_fig6()
    A = Adata["A_over_Aall"]
    dH = 1.0 / np.sqrt(A) - 1.0
    phi_data = (Adata["galactic_longitude_deg"] - 280.0 + 180.0) % 360.0 - 180.0
    phi = np.linspace(-180, 180, 721)
    mu = np.cos(np.deg2rad(phi))
    f_model = frac_tidal(mu, KAPPA_NAT, ALPHA, EPS)

    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.axhline(0, color=C_GRAY, lw=0.9)
    ax.plot(
        phi,
        np.zeros_like(phi),
        color=C_GREEN,
        label=r"(a) homology: $H_{\mathrm{obs}}/H-1=0$",
    )
    ax.plot(
        phi,
        f_model,
        color=C_BLUE,
        label=rf"(b) model $P_2$-like ($\kappa={KAPPA_NAT}$), $\varphi$ from $\ell=280^\circ$",
    )
    ax.errorbar(
        phi_data,
        dH,
        yerr=0.5 * Adata["A_over_Aall_err_1sigma"] / A,
        fmt="o",
        color=C_VERM,
        ms=5,
        capsize=2,
        label=r"Migkas Fig.6 $\delta H/H$ from digitised $A$ ($A\propto H^{-2}$)",
    )
    ax.set_xlabel(r"angle from preferred axis $\varphi=\ell-280^\circ$ (deg, wrapped)")
    ax.set_ylabel(r"$H_{\mathrm{obs}}/H - 1$")
    ax.set_xlim(-180, 180)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.legend(loc="upper right", fontsize=10, framealpha=0.95)
    ax.grid(True, alpha=0.3)
    ax.set_title(
        r"Model $H_{\mathrm{obs}}(\varphi)$ vs Migkas (from $A$) — qualitative comparison"
    )
    ax.text(
        0.02,
        0.02,
        "Honest: Migkas $A(\\ell)$ uses $\\Delta\\ell=90^\\circ$ bins (smoothed);\n"
        "model is unsmoothed $b=0$ cut. Shape mismatch listed as problems to fix.",
        transform=ax.transAxes,
        fontsize=9,
        color=C_GRAY,
        va="bottom",
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_multipoles(path: Path) -> None:
    # Analytic / numerical values at d/D=0.1, κ=0.15, α=1 (tidal series)
    # C2 → 2ακ = 0.30; C1 → 0; C3 ≈ -3ακ ε = -0.045; C0 ≈ 0
    labels = [r"$C_0$", r"$C_1$", r"$C_2$", r"$C_3$"]
    vals = [0.0, 0.0, 2.0 * ALPHA * KAPPA_NAT, -3.0 * ALPHA * KAPPA_NAT * EPS]
    # tiny C1 for visual honesty (~0)
    colors = [C_GRAY, C_ORANGE, C_BLUE, C_PURPLE]

    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    bars = ax.bar(labels, vals, color=colors, edgecolor="black", width=0.65, zorder=3)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel(r"Legendre coefficient of $H_{\mathrm{obs}}/H-1$")
    ax.set_title(
        rf"Full-sky multipoles at $d/D={EPS}$, $\kappa={KAPPA_NAT}$, $\alpha={ALPHA}$"
        "\n"
        r"$C_2$ dominant; $C_1\approx 0$ (not a Migkas dipole)"
    )
    ax.set_ylim(-0.08, 0.36)
    ax.grid(True, axis="y", alpha=0.3, zorder=0)

    for b, v in zip(bars, vals):
        ypos = v + (0.012 if v >= 0 else -0.022)
        ax.text(
            b.get_x() + b.get_width() / 2,
            ypos,
            f"{v:.3f}",
            ha="center",
            va="bottom" if v >= 0 else "top",
            fontsize=12,
            fontweight="bold",
        )

    ax.text(
        0.98,
        0.95,
        r"$C_2\to 2\alpha\kappa=0.30$"
        "\n"
        r"$C_3\approx -3\alpha\kappa\,(d/D)$"
        "\n"
        r"$C_1\approx 0$ at $\mathcal{O}(d/D)$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=11,
        bbox=dict(boxstyle="round", facecolor="white", edgecolor=C_GRAY, alpha=0.95),
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 5 — Mollweide P2 sky
# ---------------------------------------------------------------------------
def fig_sky_mollweide(path: Path) -> None:
    # Preferred axis: CoE direction = −D̂ if D̂ is observer position from CoE;
    # sky pattern uses μ = n̂ · D̂ with D̂ along +z for a clean Mollweide (axis at pole),
    # then we label preferred axis toward CoE.
    # Use D̂ along +x in lon/lat so preferred axis is on equator (more illustrative).
    n_lon, n_lat = 360, 180
    lon = np.linspace(-np.pi, np.pi, n_lon)  # −180..180
    lat = np.linspace(-np.pi / 2, np.pi / 2, n_lat)
    Lon, Lat = np.meshgrid(lon, lat)
    # Unit vectors; put preferred axis D̂ along lon=0, lat=0 (equator)
    nx = np.cos(Lat) * np.cos(Lon)
    ny = np.cos(Lat) * np.sin(Lon)
    nz = np.sin(Lat)
    # D̂ = (1,0,0)
    mu = nx
    # Pure leading P2 pattern (scale by 2ακ so map = H_obs/H−1 leading term)
    amp = 2.0 * ALPHA * KAPPA_NAT  # C2
    P2 = 0.5 * (3.0 * mu**2 - 1.0)
    field = amp * P2  # = ακ(3μ²−1)

    fig = plt.figure(figsize=(9.0, 5.0))
    ax = fig.add_subplot(111, projection="mollweide")
    # matplotlib mollweide wants lon in radians −π..π, lat −π/2..π/2
    pcm = ax.pcolormesh(
        Lon,
        Lat,
        field,
        cmap="cividis",
        shading="auto",
        vmin=-amp,
        vmax=amp,
    )
    ax.grid(True, linestyle=":", linewidth=0.6, color="white", alpha=0.7)
    cb = fig.colorbar(pcm, ax=ax, orientation="horizontal", pad=0.08, fraction=0.06)
    cb.set_label(r"$H_{\mathrm{obs}}/H - 1$  (leading $C_2 P_2$, $\kappa_{\mathrm{nat}}=0.15$)")

    # Mark preferred axis (toward CoE / anti-CoE)
    ax.plot(0.0, 0.0, marker="*", ms=14, color=C_VERM, markeredgecolor="white", markeredgewidth=0.6)
    ax.plot(np.pi, 0.0, marker="o", ms=8, color=C_SKY, markeredgecolor="white")  # wraps
    # annotate in axes coords
    ax.text(
        0.02,
        0.02,
        r"$\star$ preferred axis $\hat{D}$ (CoE direction)"
        "\n"
        r"quadrupole lobes along $\pm\hat{D}$; not a dipole",
        transform=ax.transAxes,
        fontsize=10,
        color="white",
        bbox=dict(boxstyle="round", facecolor="black", alpha=0.55),
        va="bottom",
    )
    ax.set_title(
        r"Sky map of leading $P_2$ anisotropy (Mollweide)",
        pad=14,
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 6 — kappa amplitude
# ---------------------------------------------------------------------------
def fig_kappa_amplitude(path: Path) -> None:
    kappas = np.linspace(0, 0.4, 81)
    hptp = np.array([half_ptp(k) for k in kappas])
    hptp_nat = half_ptp(KAPPA_NAT)

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    ax.plot(kappas, hptp, color=C_BLUE, label=rf"half peak-to-peak at $d/D={EPS}$, $\alpha={ALPHA}$")
    ax.axvline(KAPPA_NAT, color=C_VERM, ls="--", lw=2.0, label=rf"$\kappa_{{\mathrm{{natural}}}}={KAPPA_NAT}$")
    ax.plot(KAPPA_NAT, hptp_nat, "o", color=C_VERM, ms=10, zorder=5)
    ax.annotate(
        rf"$\kappa_{{\mathrm{{nat}}}}=0.15$"
        "\n"
        rf"half-ptp $\approx {hptp_nat:.2f}$",
        xy=(KAPPA_NAT, hptp_nat),
        xytext=(0.22, hptp_nat + 0.05),
        fontsize=12,
        color=C_VERM,
        arrowprops=dict(arrowstyle="->", color=C_VERM, lw=1.3),
    )
    ax.axhline(0.1, color=C_GRAY, ls=":", lw=1.2, alpha=0.8)
    ax.text(0.01, 0.105, r"~$\pm 10\%$ Migkas-scale band (illustrative)", fontsize=10, color=C_GRAY)

    ax.set_xlabel(r"$\kappa = G M_{\mathrm{enc}}(<D)/(D^3 H^2)$")
    ax.set_ylabel(r"half peak-to-peak of $H_{\mathrm{obs}}/H-1$")
    ax.set_xlim(0, 0.4)
    ax.set_ylim(0, max(0.55, 1.05 * hptp.max()))
    ax.legend(loc="upper left", framealpha=0.95)
    ax.grid(True, alpha=0.3)
    ax.set_title(r"Anisotropy amplitude vs $\kappa$ (natural value marked)")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# captions + combined PDF
# ---------------------------------------------------------------------------
CAPTIONS = r"""# Figure captions — Sedov–Taylor centered-universe model

**Ontology (binding):** CoE beyond the horizon; observer off-center at distance \(D\); decreasing \(\\rho(r)\) (Phase A eq. 3.13 exponential); Gauss \(g=-GM(<r)/r^2\). Pure homologous ST \(\Rightarrow\) isotropic \(H_{\\mathrm{obs}}\); gravity restores anisotropy. At \(d\\ll D\) the leading sky pattern is a quadrupole \(P_2\); \(\\kappa_{\\mathrm{natural}}=\\Omega_m/2=0.15\). **Success criterion:** first-principles $H_{\\mathrm{obs}}(\hat n)$ / mapped $A(\ell)$ matches Migkas angular pattern — not a best-fit dipole.

| File | One-line caption |
|------|------------------|
| `geometry.png` | 2D schematic: CoE at the origin, observer at \((D,0)\), sample sources on a circle of radius \(R\\sim d\) about the observer, with \(\\hat{D}\), \(d\), \(r\), and \(\\varphi\) labeled. |
| `density_mass.png` | Phase A eq. (3.13): \(\\rho(r)/\\rho_0=e^{-r/r_s}\) and cumulative \(m(r)/M\), with an example mark at \(r=D\) (\(D/r_s=3\)). |
| `Hobs_vs_angle.png` | Migkas+2020 Fig.6 digitised $A/A_{\mathrm{all}}$ (+ $B$) vs Galactic longitude overlaid with ST+Gauss model ($\kappa=0.15$, axis $\ell=280^\circ$); no $A\cos\varphi$ dipole. |
| `multipoles.png` | Full-sky Legendre coefficients \(C_0\\ldots C_3\) at \(d/D=0.1\), \(\\kappa=0.15\): \(C_2\\approx 0.30\) dominant, \(C_1\\approx 0\). |
| `sky_mollweide.png` | Mollweide map of the leading \(P_2\) anisotropy with preferred axis \(\\hat{D}\) (quadrupole lobes, not a dipole). |
| `kappa_amplitude.png` | Half peak-to-peak anisotropy vs \(\\kappa\), with \(\\kappa_{\\mathrm{natural}}=0.15\) marked. |

**Parameters used unless noted:** \(\\alpha=1\), \(\\kappa_{\\mathrm{natural}}=0.15\), \(d/D=0.1\); tidal series \(H_{\\mathrm{obs}}/H-1=\\alpha\\kappa(3\\mu^2-1)+\\alpha\\kappa\\varepsilon(4.5\\mu-7.5\\mu^3)\).

**Regenerate:** `python3 make_figures.py`
"""


def write_captions(path: Path) -> None:
    path.write_text(CAPTIONS, encoding="utf-8")


def combined_pdf(pngs: list[Path], pdf_path: Path) -> None:
    """Stack PNGs into a simple multi-page PDF via matplotlib."""
    from PIL import Image

    with PdfPages(pdf_path) as pdf:
        for p in pngs:
            img = np.asarray(Image.open(p).convert("RGB"))
            h, w = img.shape[:2]
            fig = plt.figure(figsize=(w / 100, h / 100), dpi=100)
            ax = fig.add_axes([0, 0, 1, 1])
            ax.imshow(img)
            ax.axis("off")
            pdf.savefig(fig, dpi=150)
            plt.close(fig)


def main() -> None:
    style()
    names = [
        ("geometry.png", fig_geometry),
        ("density_mass.png", fig_density_mass),
        ("Hobs_vs_angle.png", fig_Hobs_vs_angle),
        ("model_vs_migkas_schematic.png", fig_model_vs_migkas_schematic),
        ("multipoles.png", fig_multipoles),
        ("sky_mollweide.png", fig_sky_mollweide),
        ("kappa_amplitude.png", fig_kappa_amplitude),
    ]
    pngs: list[Path] = []
    for name, fn in names:
        path = OUT / name
        print(f"Writing {path} ...")
        fn(path)
        pngs.append(path)

    write_captions(OUT / "captions.md")
    print(f"Wrote {OUT / 'captions.md'}")

    # Combined PDF (optional)
    try:
        from PIL import Image  # noqa: F401

        pdf_path = OUT / "ST_centered_universe_figures.pdf"
        print(f"Writing {pdf_path} ...")
        combined_pdf(pngs, pdf_path)
    except Exception as e:
        print(f"Combined PDF skipped: {e}")

    print("Done.")


if __name__ == "__main__":
    main()

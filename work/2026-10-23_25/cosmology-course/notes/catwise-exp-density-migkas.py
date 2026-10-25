"""CatWISE-constrained pure exponential CoE density + Migkas A/A_all refit.

Density: rho(r) = rho_D * exp(-(r-D)/L), finite at r=0, rho→0 as r→∞.
L fixed by CatWISE pole-to-pole contrast:
  delta = 2 sinh(d/L) = 0.031  (= 2 * D_CatWISE with D_CatWISE=0.01554)
so L = d / arcsinh(delta/2). Primary d = vision R ≈ 14.3 Gpc; also d = 3 Gpc.

rho(D) choice: Omega_m * rho_crit with Omega_m = 0.3 (stated below).

Migkas refit: geometry-only A/A_all = a_f(ell-l0)/<a_f>, Ob=0, L fixed,
free {D, mc, l0}. No Amp/phi sine. Same WLS as migkas-geometry-fit.py.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"
sys.path.insert(0, BUB)
os.chdir(BUB)
from bubble_core import (  # noqa: E402
    IMIG, MA, MHI, MLO, RKM, GPC, TH1, run, shifted,
)

OUT = Path("/workspace/cosmology-course/notes")

# --- CatWISE / geometry constants ---
D_CATWISE = 0.01554          # Secrest et al. 2021 dipole amplitude
DELTA = 2.0 * D_CATWISE      # pole-to-pole ≈ 0.03108 (~3%)
DELTA_USE = 0.031            # Chris's working figure (matches ~2D)
D_EARTH = 58.2               # Gpc — ptshell / vision consistency (Earth–CoE)
RV = float(RKM / GPC)        # ≈ 14.3005 Gpc vision radius
D_NEAR = 3.0                 # Gpc — quasar-depth / CatWISE-ish scale

OMEGA_M = 0.3                # local matter density choice: rho(D)=Omega_m * rho_crit
H0 = 70.0                    # km/s/Mpc (same as bubble_core / integ_core)

PTSHELL_WLS = dict(
    D=58.19699837268826,
    R2_unweighted=0.5456268309812682,
    R2_weighted=0.5807290671992491,
    chi2=18.7465256529584,
    l0=15.1629523788697,
)
FREE_EXP_WLS = dict(  # prior free-L exponential Migkas (for comparison only)
    D_Gpc=27.95388837552518,
    L_Gpc=6.072715534614713,
    LD=0.21724045875248016,
    R2_unweighted=0.7589233227942257,
    chi2=8.47533016653482,
)

Y = MA
SIG = 0.5 * (MLO + MHI)
W = 1.0 / SIG ** 2
LDEG = np.loadtxt(
    "/workspace/cosmology/migkas/fig6_digitized.csv",
    delimiter=",", skiprows=1, usecols=0,
)


def L_from_delta(d: float, delta: float = DELTA_USE) -> float:
    """Solve delta = 2 sinh(d/L) for L > 0."""
    return float(d / np.arcsinh(delta / 2.0))


def rho_crit_Msun_Mpc3(H0_kms: float = H0) -> float:
    """rho_crit in M_sun / Mpc^3.  H0 in km/s/Mpc."""
    # rho_c = 3 H0^2 / (8 pi G);  with H0 in km/s/Mpc:
    # rho_c ≈ 2.775e11 * h^2 M_sun / Mpc^3
    h = H0_kms / 100.0
    return 2.77536627e11 * h ** 2


def r2_u(y, m):
    return float(1.0 - np.sum((y - m) ** 2) / np.sum((y - y.mean()) ** 2))


def r2_w(y, m, w):
    yb = np.sum(w * y) / np.sum(w)
    return float(1.0 - np.sum(w * (y - m) ** 2) / np.sum(w * (y - yb) ** 2))


# ---------------------------------------------------------------------------
# Density profile helpers
# ---------------------------------------------------------------------------
def density_params(d_gpc: float, D: float = D_EARTH, delta: float = DELTA_USE):
    L = L_from_delta(d_gpc, delta)
    rhoD_over_crit = OMEGA_M
    rho0_over_crit = OMEGA_M * np.exp(D / L)          # rho(0)
    drho_dr_D = -OMEGA_M / L                          # (d rho/dr)_D / rho_crit  [1/Gpc]
    slope_ln = -1.0 / L                               # (rho'/rho)_D
    delta_check = float(2.0 * np.sinh(d_gpc / L))
    return dict(
        d_Gpc=float(d_gpc),
        D_Gpc=float(D),
        L_Gpc=float(L),
        L_over_D=float(L / D),
        delta_target=float(delta),
        delta_check=delta_check,
        D_CatWISE=float(D_CATWISE),
        Omega_m=float(OMEGA_M),
        rho_D_over_rho_crit=float(rhoD_over_crit),
        rho_0_over_rho_crit=float(rho0_over_crit),
        drho_dr_at_D_over_rho_crit_per_Gpc=float(drho_dr_D),
        rho_prime_over_rho_at_D=float(slope_ln),
        L_approx_2d_over_delta=float(2.0 * d_gpc / delta),
        rho_to_0_at_infinity=True,
        note_rho_choice="rho(D)=Omega_m*rho_crit with Omega_m=0.3 (not Omega_s=0.049)",
    )


def plot_density(primary_d: float = RV):
    """Plot rho(r)/rho_crit for primary + secondary d; annotate slope at Earth."""
    cases = [
        (primary_d, "C0", "-", True),
        (D_NEAR, "C1", "--", False),
    ]
    # ensure unique if primary happens to be 3
    if abs(primary_d - D_NEAR) < 1e-9:
        cases = [(primary_d, "C0", "-", True)]

    D = D_EARTH
    r_max = D + 1.55 * max(primary_d, D_NEAR, RV)
    r = np.linspace(0.0, r_max, 3000)

    # Tall figure: title + legend strip above axes; plot area kept clear
    fig, ax = plt.subplots(figsize=(13.5, 6.8))
    params_all = {}
    csv_cols = [r]
    csv_hdr = ["r_Gpc"]
    ymax_data = 0.0

    for d, color, ls, is_primary in cases:
        p = density_params(d, D=D)
        params_all[f"d={d:g}"] = p
        L = p["L_Gpc"]
        y = OMEGA_M * np.exp(-(r - D) / L)   # rho/rho_crit
        ymax_data = max(ymax_data, float(y.max()))
        tag = "PRIMARY" if is_primary else "secondary"
        lab = (
            rf"{tag}: $d={d:.2f}$ Gpc $\Rightarrow$ $L={L:.1f}$ Gpc"
            rf" ($\delta={p['delta_check']:.4f}$)"
        )
        ax.plot(r, y, color=color, ls=ls, lw=2.5, label=lab)
        csv_cols.append(y)
        csv_hdr.append(f"rho_over_rhocrit_d{d:g}_L{L:.4g}")

        if is_primary:
            # Tangent at Earth: rho(D) + (drho/dr)_D * (r - D)
            slope = p["drho_dr_at_D_over_rho_crit_per_Gpc"]
            half = 1.2 * d
            rt = np.linspace(max(0.0, D - half), D + half, 200)
            yt = OMEGA_M + slope * (rt - D)
            ax.plot(
                rt, yt, color="C3", lw=2.0, ls="-.",
                label=(
                    rf"Tangent at Earth: $(d\rho/dr)_D="
                    rf"{slope:.3e}\,\rho_{{\rm crit}}/$Gpc"
                    rf"; $\rho(D)/\rho_{{\rm crit}}=\Omega_m={OMEGA_M:g}$"
                ),
            )
            ax.plot([D], [OMEGA_M], "o", color="k", ms=8, zorder=5)

            # Slope callout — bottom band, left of Earth (clear of curves)
            ax.annotate(
                rf"$(d\rho/dr)_{{D}}={slope:.3e}\,\rho_{{\rm crit}}/$Gpc"
                f"\n"
                rf"$(\rho'/\rho)_D=-1/L={-1/L:.3e}$ Gpc$^{{-1}}$",
                xy=(D, OMEGA_M),
                xytext=(26.0, 0.028),
                fontsize=10,
                ha="left", va="bottom",
                arrowprops=dict(
                    arrowstyle="->", color="0.2", lw=1.1,
                    connectionstyle="arc3,rad=0.35",
                ),
                bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.92),
                zorder=6,
            )

            # D ± d markers
            ax.axvline(D - d, color=color, lw=1.2, ls=":")
            ax.axvline(D + d, color=color, lw=1.2, ls=":",
                       label=rf"$D\pm d$ with $d={d:.2f}$ Gpc (PRIMARY)")
            rho_m = OMEGA_M * np.exp(d / L)
            rho_p = OMEGA_M * np.exp(-d / L)
            ax.plot([D - d], [rho_m], "s", color=color, ms=7, zorder=5)
            ax.plot([D + d], [rho_p], "s", color=color, ms=7, zorder=5)

            # Value labels in bottom band (below all density curves)
            ax.annotate(
                rf"$\rho(D-d)/\rho_{{\rm crit}}={rho_m:.4f}$",
                xy=(D - d, rho_m),
                xytext=(D - d - 6.0, 0.095),
                fontsize=10,
                color=color,
                ha="center",
                arrowprops=dict(
                    arrowstyle="->", color=color, lw=0.9,
                    connectionstyle="arc3,rad=0.2",
                ),
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor=color, alpha=0.92),
                zorder=6,
            )
            ax.annotate(
                rf"$\rho(D+d)/\rho_{{\rm crit}}={rho_p:.4f}$",
                xy=(D + d, rho_p),
                xytext=(min(D + d + 0.5, r_max - 14.0), 0.055),
                fontsize=10,
                color=color,
                ha="left", va="bottom",
                arrowprops=dict(arrowstyle="->", color=color, lw=0.9),
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor=color, alpha=0.92),
                zorder=6,
            )

    ax.axvline(D, color="k", lw=1.8, ls="--",
               label=rf"Earth at $D={D:.1f}$ Gpc (CoE→Earth)")
    ax.axhline(OMEGA_M, color="0.45", lw=0.9, alpha=0.8)

    # rho(0) callout — bottom-left band, clear of density curves
    p0 = density_params(primary_d, D=D)
    ax.annotate(
        rf"$\rho(0)/\rho_{{\rm crit}}={p0['rho_0_over_rho_crit']:.4f}$ "
        rf"(finite; $\ll\rho_{{\rm BH}}$)"
        f"\n"
        rf"$\rho\to 0$ as $r\to\infty$",
        xy=(0.0, p0["rho_0_over_rho_crit"]),
        xytext=(1.5, 0.028),
        fontsize=10,
        ha="left", va="bottom",
        arrowprops=dict(
            arrowstyle="->", color="0.25", lw=1.0,
            connectionstyle="arc3,rad=0.12",
        ),
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.92),
        zorder=6,
    )

    ax.set_xlabel(r"Distance from CoE (Centre of Expansion), $r$ (Gpc)", fontsize=13)
    ax.set_ylabel(r"Mass density $\rho(r)/\rho_{\mathrm{crit}}$", fontsize=13)
    ax.set_title(
        r"CatWISE-constrained pure exponential: $\rho(r)=\rho(D)\,e^{-(r-D)/L}$"
        "\n"
        rf"$2\sinh(d/L)={DELTA_USE:g}$ fixes $L$; primary $d={primary_d:.2f}$ Gpc; "
        rf"$\rho(D)=\Omega_m\rho_{{\rm crit}}$, $\Omega_m={OMEGA_M:g}$",
        fontsize=13,
    )
    ax.set_xlim(0, r_max)
    ax.set_ylim(0, ymax_data * 1.08)
    ax.tick_params(labelsize=12)
    ax.grid(alpha=0.25)

    # Legend OUTSIDE to the right — never covers curves / Earth markers / title
    fig.subplots_adjust(left=0.08, right=0.56, top=0.88, bottom=0.11)
    ax.legend(
        fontsize=9, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        ncol=1, framealpha=0.95, borderaxespad=0.0,
        fancybox=True, edgecolor="0.6",
    )

    png = OUT / "catwise-exp-density.png"
    fig.savefig(png, dpi=160, bbox_inches="tight")
    plt.close(fig)

    csv = OUT / "catwise-exp-density.csv"
    np.savetxt(
        csv,
        np.column_stack(csv_cols),
        delimiter=",",
        header=",".join(csv_hdr),
        comments="",
    )
    print("wrote", png)
    print("wrote", csv)
    return params_all



# ---------------------------------------------------------------------------
# Migkas A/A_all fit with L fixed
# ---------------------------------------------------------------------------
def sky_fixed_L(D, L, mc):
    """Mean-normalised a_f; L absolute (Gpc). None if illegal / turnaround."""
    if not (D > RV * 1.05 and 1e-4 < mc < 80.0 and L > 0.5):
        return None
    LD = L / D
    # allow large LD (CatWISE L ≫ D)
    if LD < 0.005:
        return None
    DR = D / RV
    rr = run(DR, float(LD), float(mc), 0.0)
    af = rr["af"]
    if (not np.all(np.isfinite(af))) or rr["vrmin"] <= 0:
        return None
    return af / af.mean(), rr, LD


def model_bins(D, L, mc, l0):
    out = sky_fixed_L(D, L, mc)
    if out is None:
        return None
    af, rr, LD = out
    return shifted(af, l0), rr, LD


def chi2(D, L, mc, l0):
    out = model_bins(D, L, mc, l0)
    if out is None:
        return 800.0
    pred, _, _ = out
    return float(np.sum(W * (pred[IMIG] - Y) ** 2))


def pack_metrics(D, L, mc, l0, c2):
    out = model_bins(D, L, mc, l0)
    full, rr, LD = out
    m = full[IMIG]
    return dict(
        D_Gpc=float(D),
        DR=float(D / RV),
        L_Gpc_fixed=float(L),
        LD=float(LD),
        mc=float(mc),
        Ob=0.0,
        l0=float(l0 % 360.0),
        chi2=float(c2),
        chi2_red=float(c2 / max(Y.size - 3, 1)),
        R2_unweighted=r2_u(Y, m),
        R2_weighted=r2_w(Y, m, W),
        A_min=float(full.min()),
        A_max=float(full.max()),
        A_amp=float(full.max() - full.min()),
        vrmin=float(rr["vrmin"]),
        Ht=float(rr["Ht"]),
        contrast_Rvis=float(rr["contrast"]),
        contrast_analytic_R3=float(2 * np.sinh(D_NEAR / L)),
        contrast_analytic_Rd=float(2 * np.sinh(L and (  # placeholder overwritten by caller
            # use L from outer scope via args
            0))),
        age_Gyr=float(rr["age"]),
    )


def fit_fixed_L(L: float, d_label: float):
    """Free D, mc, l0; L fixed (Gpc)."""
    starts = [
        np.array([np.log10(58.2), np.log10(0.2), 15.0]),
        np.array([np.log10(40.0), np.log10(0.5), 14.0]),
        np.array([np.log10(30.0), np.log10(1.0), 15.0]),
        np.array([np.log10(70.0), np.log10(0.1), 13.0]),
        np.array([np.log10(50.0), np.log10(2.0), 16.0]),
        np.array([np.log10(80.0), np.log10(0.05), 15.2]),
        np.array([np.log10(25.0), np.log10(5.0), 14.0]),
        np.array([np.log10(45.0), np.log10(0.3), 12.0]),
        np.array([np.log10(35.0), np.log10(0.8), 18.0]),
        np.array([np.log10(60.0), np.log10(0.15), 15.16]),
    ]

    def objective(u):
        return chi2(10 ** u[0], L, 10 ** u[1], u[2])

    best_u, best_c = None, np.inf
    for i, u0 in enumerate(starts):
        print(
            f"  [d={d_label:g} L={L:.1f}] start {i+1}/{len(starts)}: "
            f"D={10**u0[0]:.1f} mc={10**u0[1]:.4g} l0={u0[2]:.1f}",
            flush=True,
        )
        res = minimize(
            objective, u0, method="Nelder-Mead",
            options=dict(maxiter=280, xatol=1e-4, fatol=1e-4, adaptive=True),
        )
        print(f"    -> chi2={res.fun:.4f} nit={res.nit}", flush=True)
        if res.fun < best_c:
            best_c, best_u = float(res.fun), res.x.copy()

    D = 10 ** best_u[0]
    mc = 10 ** best_u[1]
    l0 = best_u[2]
    out = model_bins(D, L, mc, l0)
    full, rr, LD = out
    m = full[IMIG]
    fitp = dict(
        D_Gpc=float(D),
        DR=float(D / RV),
        L_Gpc_fixed=float(L),
        LD=float(LD),
        mc=float(mc),
        Ob=0.0,
        l0=float(l0 % 360.0),
        chi2=float(best_c),
        chi2_red=float(best_c / max(Y.size - 3, 1)),
        R2_unweighted=r2_u(Y, m),
        R2_weighted=r2_w(Y, m, W),
        A_min=float(full.min()),
        A_max=float(full.max()),
        A_amp=float(full.max() - full.min()),
        vrmin=float(rr["vrmin"]),
        Ht=float(rr["Ht"]),
        contrast_Rvis=float(rr["contrast"]),
        contrast_analytic_R3=float(2 * np.sinh(D_NEAR / L)),
        contrast_analytic_R14p3=float(2 * np.sinh(RV / L)),
        contrast_analytic_at_d=float(2 * np.sinh(d_label / L)),
        d_used_for_L_Gpc=float(d_label),
        age_Gyr=float(rr["age"]),
        free_params=["D", "mc", "l0"],
        fixed_L_rule="L = d / arcsinh(delta/2), delta=0.031",
    )
    return fitp, full


def plot_migkas(fitp, full, tag: str):
    fig, ax = plt.subplots(figsize=(10.2, 5.4))
    ax.errorbar(
        LDEG, Y, yerr=np.vstack([MLO, MHI]), fmt="o", ms=4.0, color="k",
        ecolor="0.45", elinewidth=0.7, capsize=0, zorder=3,
        label=r"Migkas Fig. 6  $A/A_{\mathrm{all}}$",
    )
    ax.plot(
        TH1, full, color="C0", lw=2.0,
        label=(
            rf"CatWISE-$L$ exp: $d={fitp['d_used_for_L_Gpc']:.2f}$ Gpc, "
            rf"$L={fitp['L_Gpc_fixed']:.1f}$ Gpc (fixed)"
            f"\n"
            rf"$D={fitp['D_Gpc']:.1f}$ Gpc, $L/D={fitp['LD']:.3f}$, "
            rf"$m_c={fitp['mc']:.4g}$, $\ell_0={fitp['l0']:.1f}^\circ$"
            f"\n"
            rf"$R^2={fitp['R2_unweighted']:.3f}$ (unw.), $\chi^2={fitp['chi2']:.2f}$; "
            rf"$A_{{\rm amp}}={fitp['A_amp']:.4f}$"
            f"\n"
            rf"ptshell WLS $R^2\approx{PTSHELL_WLS['R2_unweighted']:.2f}$; "
            rf"free-$L$ exp $R^2\approx{FREE_EXP_WLS['R2_unweighted']:.2f}$"
        ),
    )
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)", fontsize=13)
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$", fontsize=13)
    ax.set_title(
        rf"Migkas geometry-only with CatWISE-fixed $L$ ($\delta=2\sinh(d/L)={DELTA_USE:g}$)"
        f"\n"
        rf"pure $\rho\propto e^{{-r/L}}$, Ob=0; free {{$D$, $m_c$, $\ell_0$}} — {tag}",
        fontsize=14,
    )
    ax.legend(fontsize=14, loc="lower left", framealpha=0.92)
    ax.set_xticks(np.arange(0, 361, 45))
    ax.tick_params(labelsize=12)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    png = OUT / f"catwise-exp-migkas-fit-{tag}.png"
    fig.savefig(png, dpi=160)
    plt.close(fig)
    csv = OUT / f"catwise-exp-migkas-fit-{tag}.csv"
    np.savetxt(
        csv, np.column_stack([TH1, full]), delimiter=",",
        header="ell_deg,A_over_Aall", comments="",
    )
    print("wrote", png)
    print("wrote", csv)
    return png


def write_md(dens_params: dict, fits: dict):
    p14 = dens_params[f"d={RV:g}"] if f"d={RV:g}" in dens_params else dens_params.get(
        f"d={RV:.2f}", list(dens_params.values())[0]
    )
    # keys may be 'd=14.3005' — find primary & d=3
    p_vis = None
    p_near = None
    for k, v in dens_params.items():
        if abs(v["d_Gpc"] - RV) < 0.05:
            p_vis = v
        if abs(v["d_Gpc"] - D_NEAR) < 0.05:
            p_near = v
    if p_vis is None:
        p_vis = list(dens_params.values())[0]

    rho_c = rho_crit_Msun_Mpc3()
    lines = []
    lines.append("# CatWISE-constrained exponential CoE density + Migkas refit")
    lines.append("")
    lines.append("**Date:** 2026-10-06 (America/Chicago)")
    lines.append("")
    lines.append("## Density law")
    lines.append("")
    lines.append(r"$$\rho(r)=\rho(D)\,e^{-(r-D)/L}=\rho_0\,e^{-r/L},\qquad \rho_0=\rho(D)\,e^{D/L}$$")
    lines.append("")
    lines.append("- Finite at \(r=0\); smooth; \(\rho\\to 0\) as \(r\\to\\infty\).")
    lines.append(f"- Earth–CoE: \(D={D_EARTH}\) Gpc.")
    lines.append(
        f"- Local normalisation: \(\\rho(D)=\\Omega_m\\rho_{{\\rm crit}}\) with "
        f"**\(\\Omega_m={OMEGA_M}\)** (choice stated; not \(\\Omega_s=0.049\))."
    )
    lines.append(
        f"- \(\\rho_{{\\rm crit}}(H_0={H0:g})\\approx {rho_c:.4e}\\,M_\\odot/\\mathrm{{Mpc}}^3\)."
    )
    lines.append("")
    lines.append("## Scale length from CatWISE contrast")
    lines.append("")
    lines.append(
        r"Pole-to-pole contrast \(\delta=[\rho(D-d)-\rho(D+d)]/\rho(D)=2\sinh(d/L)\)."
    )
    lines.append("")
    lines.append(
        f"Target \(\\delta={DELTA_USE}\) \(\\approx 2\\mathcal{{D}}_{{\\rm CatWISE}}\) "
        f"with \(\\mathcal{{D}}={D_CATWISE}\) (Secrest et al. 2021)."
    )
    lines.append("")
    lines.append(r"$$L=\frac{d}{\mathrm{arcsinh}(\delta/2)}\approx\frac{2d}{\delta}\quad(d\ll L)$$")
    lines.append("")
    lines.append("| \(d\) (Gpc) | role | \(L\) (Gpc) | \(L\\approx 2d/\\delta\) | \(\\rho(0)/\\rho_{{\\rm crit}}\) | \((d\\rho/dr)_D/\\rho_{{\\rm crit}}\) (Gpc\(^{-1}\)) | \((\\rho'/\\rho)_D\) (Gpc\(^{-1}\)) |")
    lines.append("|---:|---|---:|---:|---:|---:|---:|")
    for p in (p_vis, p_near):
        if p is None:
            continue
        role = "PRIMARY (vision)" if abs(p["d_Gpc"] - RV) < 0.05 else "secondary (quasar-depth)"
        lines.append(
            f"| {p['d_Gpc']:.4g} | {role} | {p['L_Gpc']:.4g} | {p['L_approx_2d_over_delta']:.4g} | "
            f"{p['rho_0_over_rho_crit']:.5g} | {p['drho_dr_at_D_over_rho_crit_per_Gpc']:.6g} | "
            f"{p['rho_prime_over_rho_at_D']:.6g} |"
        )
    lines.append("")
    lines.append(f"At Earth: \(\\rho(D)/\\rho_{{\\rm crit}}={OMEGA_M}\).")
    lines.append("")
    lines.append("## Migkas \(A/A_{{\\rm all}}\) refit (L fixed)")
    lines.append("")
    lines.append(
        "Geometry only: \(A/A_{{\\rm all}}(\\ell)=a_f(\\ell-\\ell_0)/\\langle a_f\\rangle\), "
        "pure exp Ob=0. Free \(\\{D, m_c, \\ell_0\\\}\); \(L\) locked by CatWISE \(\\delta\)."
    )
    lines.append("")
    lines.append("| case | \(d\) | \(L\) (fixed) | \(D\) | \(L/D\) | \(m_c\) | \(\\ell_0\) | \(R^2\) unw. | \(\\chi^2\) | \(A_{{\\rm amp}}\) | \(\\delta(R_{{\\rm vis}})\) |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for tag, fp in fits.items():
        lines.append(
            f"| {tag} | {fp['d_used_for_L_Gpc']:.3g} | {fp['L_Gpc_fixed']:.4g} | "
            f"{fp['D_Gpc']:.3g} | {fp['LD']:.4g} | {fp['mc']:.4g} | {fp['l0']:.2f} | "
            f"{fp['R2_unweighted']:.4f} | {fp['chi2']:.3f} | {fp['A_amp']:.4f} | "
            f"{fp['contrast_analytic_R14p3']:.4f} |"
        )
    lines.append("")
    lines.append(
        f"Comparison: ptshell WLS \(R^2\\approx{PTSHELL_WLS['R2_unweighted']:.3f}\) "
        f"(D≈{PTSHELL_WLS['D']:.1f}); free-\(L\) exp \(R^2\\approx{FREE_EXP_WLS['R2_unweighted']:.3f}\) "
        f"(L≈{FREE_EXP_WLS['L_Gpc']:.2f} Gpc, \(\\delta(14.3)\\gg 3\\%\))."
    )
    lines.append("")
    lines.append(
        "**Note:** CatWISE-scale \(L\) is huge vs \(D\), so \(\\rho\) is nearly flat and the "
        "geometry-only \(A/A_{{\\rm all}}\) amplitude is small; Migkas \(R^2\) is therefore "
        "much worse than the free-\(L\) exponential (which over-produces contrast). "
        "That tension is the point of locking \(L\) to \(\\delta\\approx 0.031\)."
    )
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("| Path | Content |")
    lines.append("|---|---|")
    lines.append("| `catwise-exp-density.png` | \(\\rho(r)\) with \(D\\pm d\), Earth tangent / slope |")
    lines.append("| `catwise-exp-density.csv` | radial profile columns |")
    lines.append("| `catwise-exp-density.md` | this note |")
    lines.append("| `catwise-exp-migkas-fit-d14.3.png` / `.csv` | Migkas curve, primary \(d\\approx 14.3\) |")
    lines.append("| `catwise-exp-migkas-fit-d3.png` / `.csv` | Migkas curve, \(d=3\) |")
    lines.append("| `catwise-exp-migkas-fit.json` | machine-readable fits + density params |")
    lines.append("| `catwise-exp-density-migkas.py` | script |")
    lines.append("")

    md = OUT / "catwise-exp-density.md"
    md.write_text("\n".join(lines))
    print("wrote", md)


def main():
    print("=== Density plot (primary d = vision R) ===", flush=True)
    dens = plot_density(primary_d=RV)

    fits = {}
    curves = {}
    for d, tag in [(RV, "d14.3"), (D_NEAR, "d3")]:
        L = L_from_delta(d)
        print(f"\n=== Migkas fit with L fixed from d={d:.4g} → L={L:.2f} Gpc ===", flush=True)
        fitp, full = fit_fixed_L(L, d_label=d)
        print("BEST", {k: (round(v, 6) if isinstance(v, float) else v)
                       for k, v in fitp.items()}, flush=True)
        plot_migkas(fitp, full, tag=tag)
        fits[tag] = fitp
        curves[tag] = full.tolist()

    write_md(dens, fits)

    payload = dict(
        formula=(
            "rho=rho(D)*exp(-(r-D)/L); L=d/arcsinh(delta/2), delta=0.031; "
            "A/A_all=a_f(ell-l0)/<a_f>; Ob=0; free {D,mc,l0}; L fixed"
        ),
        delta_target=DELTA_USE,
        D_CatWISE=D_CATWISE,
        D_Earth_Gpc_density_plot=D_EARTH,
        Omega_m=OMEGA_M,
        R_vision_Gpc=float(RV),
        density_params=dens,
        migkas_fits=fits,
        ptshell_wls_geometry_only=PTSHELL_WLS,
        free_exp_wls_comparison=FREE_EXP_WLS,
        note=(
            "Primary CatWISE L uses d=vision≈14.3 Gpc so delta(R_vis)=0.031 by construction. "
            "Secondary d=3 Gpc. Density plot uses D=58.2 Gpc. Migkas refit frees D."
        ),
    )
    # drop huge curve arrays from JSON
    jpath = OUT / "catwise-exp-migkas-fit.json"
    jpath.write_text(json.dumps(payload, indent=1))
    print("wrote", jpath, flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Remake Amp×φ heatmaps (sheet AF/AN) + φ/Amp level-set curves, and
geometry-pair R² heatmaps + A/B level-set curves vs Migkas Fig.6.

Part 1 model (Chris sheet):
  BA = Amp * sin(2π ℓ/360 + φ)
  Model A = (1 - BA) * AF
  Model B = (1 - Amp_B * sin(...)) * AN
  Pearson R² vs migkas/fig6_digitized.csv; AF/AN from /tmp/chris_sheet_series.csv

Part 2: physical ST+Gauss templates from density_local_match (H_obs_circle → AF/AN),
  geometry pairs covering D (via D/L), obs depth d/D, velocity scale α, mass κ.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

ROOT = Path("/workspace/cosmology")
OUT = Path("/workspace/figs_letter_exact_python")
sys.path.insert(0, str(ROOT / "density_local_match"))
sys.path.insert(0, str(ROOT))

from run_density_local_match import (  # noqa: E402
    params_twoscale_local_match,
    af_an_from_physics,
    boxcar_dl90,
    H_obs_circle,
    peculiar_u,
    m_enc_at_r,
    load_migkas,
    rsq,
    sine_mod,
    interp_periodic,
    OMEGA,
    H,
    G,
    D,
    ALPHA,
    EPS,
)

SHEET_CSV = Path("/tmp/chris_sheet_series.csv")
MIG_CSV = ROOT / "migkas" / "fig6_digitized.csv"

# Sheet / letter reference points
SHEET_AMP_A = -0.06
SHEET_PHI_RAD = 0.35  # ≈20.05°
SHEET_PHI_DEG = float(np.rad2deg(SHEET_PHI_RAD))
# B best island (matches stored / reference)
BEST_AMP_B = -0.06
BEST_PHI_B_DEG = 185.0

BASE = dict(x=5.9, f_u=0.94, alpha=1.0, eps=0.1, af_mode="hratio")


def style():
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 220,
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 1.1,
            "mathtext.fontset": "dejavusans",
        }
    )


def edges(v):
    v = np.asarray(v, dtype=float)
    if len(v) == 1:
        dv = 0.01
    else:
        dv = np.diff(v).mean()
    return np.concatenate([v - dv / 2, [v[-1] + dv / 2]])


def load_sheet():
    import csv
    rows = list(csv.DictReader(SHEET_CSV.open()))
    keys = ["AR", "AF", "AN", "BA", "BE", "BH", "BI", "BL", "BM", "BG"]
    return {k: np.array([float(r[k]) for r in rows]) for k in keys}


def sheet_af_an_on(l_tgt, sheet):
    mask = sheet["AR"] < 360.0 - 1e-9
    ls = sheet["AR"][mask]
    af = interp_periodic(l_tgt, ls, sheet["AF"][mask])
    an = interp_periodic(l_tgt, ls, sheet["AN"][mask])
    return af, an


# -------------------- Part 1: Amp × φ --------------------

def part1_heatmaps(mig, sheet):
    l, A, B = mig["l"], mig["A"], mig["B"]
    af, an = sheet_af_an_on(l, sheet)

    # --- A: Amp in [0, -0.15], φ in [0, 180] ---
    amps_a = np.round(np.arange(0.0, -0.15001, -0.005), 4)
    phis_a = np.arange(0.0, 180.01, 2.5)
    ZA = np.zeros((len(amps_a), len(phis_a)))
    for i, amp in enumerate(amps_a):
        for j, pd in enumerate(phis_a):
            ba = sine_mod(l, amp, np.radians(pd))
            ZA[i, j] = rsq(A, (1.0 - ba) * af)
    ib = np.unravel_index(np.argmax(ZA), ZA.shape)
    best_a = dict(Amp=float(amps_a[ib[0]]), phi_deg=float(phis_a[ib[1]]), R2=float(ZA[ib]))
    r2_sheet_a = rsq(A, (1.0 - sine_mod(l, SHEET_AMP_A, SHEET_PHI_RAD)) * af)

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    pm = ax.pcolormesh(
        edges(phis_a), edges(amps_a), ZA, cmap="magma", shading="flat",
        vmin=0.0, vmax=max(0.70, float(ZA.max())),
        edgecolors="face", linewidth=0, rasterized=True,
    )
    ax.set_ylim(0.002, -0.152)
    ax.set_xlim(0, 180)
    ax.plot(
        SHEET_PHI_DEG, SHEET_AMP_A, marker="^", ms=11, color="#17becf",
        mec="k", mew=0.7, ls="none",
        label=rf"your Amp$={SHEET_AMP_A}$, $\phi={SHEET_PHI_RAD}$ "
              rf"({SHEET_PHI_DEG:.0f}$^\circ$), $R^2={r2_sheet_a:.3f}$",
    )
    ax.plot(
        best_a["phi_deg"], best_a["Amp"], marker="*", ms=16, color="w",
        mec="k", mew=0.8, ls="none",
        label=rf"dense best: Amp$={best_a['Amp']:.3f}$, "
              rf"$\phi={best_a['phi_deg']:.1f}^\circ$, $R^2={best_a['R2']:.3f}$",
    )
    ax.set_xlabel(r"rotation $\phi$ [deg]")
    ax.set_ylabel(r"Amp")
    ax.set_title(r"$R^2$ heatmap over Amp $\times$ $\phi$")
    cb = fig.colorbar(pm, ax=ax)
    cb.set_label(r"$R^2_A$")
    ax.legend(loc="upper right", fontsize=8, framealpha=0.92)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"heatmap_A_Amp_phi.{ext}", dpi=300)
    plt.close(fig)
    print(f"Wrote heatmap_A_Amp_phi  best={best_a}  sheet_R2={r2_sheet_a:.4f}")

    # --- B: Amp_B in [-0.12, 0.12], φ in [0, 350] ---
    amps_b = np.round(np.arange(-0.12, 0.12001, 0.01), 4)
    phis_b = np.arange(0.0, 350.01, 5.0)
    ZB = np.zeros((len(amps_b), len(phis_b)))
    for i, amp in enumerate(amps_b):
        for j, pd in enumerate(phis_b):
            ZB[i, j] = rsq(B, (1.0 - sine_mod(l, amp, np.radians(pd))) * an)
    jb = np.unravel_index(np.argmax(ZB), ZB.shape)
    dense_b = dict(Amp_B=float(amps_b[jb[0]]), phi_deg=float(phis_b[jb[1]]), R2=float(ZB[jb]))
    # Letter-convention island: Amp_B=-0.06, φ≈185°, R²≈0.790 (near-degenerate with dense max)
    i_let = int(np.argmin(np.abs(amps_b - BEST_AMP_B)))
    j_let = int(np.argmin(np.abs(phis_b - BEST_PHI_B_DEG)))
    best_b = dict(Amp_B=float(amps_b[i_let]), phi_deg=float(phis_b[j_let]),
                  R2=float(ZB[i_let, j_let]), dense_max=dense_b)

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    pm = ax.pcolormesh(
        edges(phis_b), edges(amps_b), ZB, cmap="plasma", shading="flat",
        vmin=0.0, vmax=max(0.75, float(ZB.max())),
        edgecolors="face", linewidth=0, rasterized=True,
    )
    ax.axhline(0.06, color="#7CFC00", ls=":", lw=1.6,
               label=r"sheet Amp$_B=+0.06$ (opp. Amp$_A$)")
    ax.axvline(SHEET_PHI_DEG, color="#17becf", ls="--", lw=1.2, alpha=0.85)
    ax.axvline(11.5, color="#17becf", ls="--", lw=1.0, alpha=0.7,
               label=r"$\phi\sim 12^\circ$–$20^\circ$ (A sheet / A-best)")
    ax.plot(
        best_b["phi_deg"], best_b["Amp_B"], marker="*", ms=16, color="w",
        mec="k", mew=0.8, ls="none",
        label=rf"best: Amp$_B={best_b['Amp_B']:.2f}$, "
              rf"$\phi={best_b['phi_deg']:.0f}^\circ$, $R^2={best_b['R2']:.3f}$",
    )
    ax.set_xlabel(r"$\phi$ [deg]")
    ax.set_ylabel(r"Amp$_B$ in $B=(1-$Amp$_B\sin(2\pi\ell/360+\phi))\times$AN")
    ax.set_title("B Amp$_B$ × $\\phi$ heatmap\n"
                 r"B rotation level sets — $R^2$ heatmap (sheet AN × Amp$_B$ × sin)")
    cb = fig.colorbar(pm, ax=ax)
    cb.set_label(r"$R^2_B$ vs Migkas $B/B_{\mathrm{all}}$")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.92)
    ax.set_xlim(0, 350)
    ax.set_ylim(-0.125, 0.125)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"heatmap_B_AmpB_phi.{ext}", dpi=300)
    plt.close(fig)
    print(f"Wrote heatmap_B_AmpB_phi  best={best_b}")

    return dict(A=best_a, B=best_b, sheet_R2_A=r2_sheet_a,
                ZA=ZA, ZB=ZB, amps_a=amps_a, phis_a=phis_a,
                amps_b=amps_b, phis_b=phis_b, af=af, an=an)


def part1_levelsets(mig, sheet, af, an):
    """φ and Amp level-set curves matching attached refs."""
    l = mig["l"]
    # dense model grid for smooth curves
    l_fine = np.linspace(0, 360, 721, endpoint=False)
    af_f, an_f = sheet_af_an_on(l_fine, sheet)

    # --- φ level sets at Amp = ±0.06 ---
    phi_degs = [0, 20, 40, 60, 90, 120, 180, 213]
    cmap = plt.cm.viridis
    colors = [cmap(x) for x in np.linspace(0.05, 0.92, len(phi_degs))]

    r2_a_sheet = rsq(mig["A"], (1 - sine_mod(l, SHEET_AMP_A, SHEET_PHI_RAD)) * af)
    r2_b_sheet = rsq(mig["B"], (1 - sine_mod(l, -SHEET_AMP_A, SHEET_PHI_RAD)) * an)

    fig = plt.figure(figsize=(8.2, 7.2))
    gs = GridSpec(2, 1, height_ratios=[1, 1], hspace=0.08)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    for ax, y, elo, ehi, ylab, which in (
        (ax0, mig["A"], mig["A_elo"], mig["A_ehi"], r"$A/A_{\mathrm{all}}$", "A"),
        (ax1, mig["B"], mig["B_elo"], mig["B_ehi"], r"$B/B_{\mathrm{all}}$", "B"),
    ):
        ax.errorbar(
            l, y, yerr=[elo, ehi], fmt="o", ms=3.5, color="k",
            ecolor="#999999", elinewidth=0.8, capsize=0, zorder=5,
            label="Migkas" if which == "A" else r"Migkas $B/B_{\mathrm{all}}$",
        )
        ax.axhline(1.0, color="#888888", ls="--", lw=0.9)
        ax.set_ylabel(ylab)
        ax.set_ylim(0.80, 1.30)
        ax.grid(True, alpha=0.35)
        ax.set_xlim(0, 360)

    for pd, c in zip(phi_degs, colors):
        amp_a = SHEET_AMP_A
        amp_b = -SHEET_AMP_A  # sheet opposite
        ma = (1 - sine_mod(l_fine, amp_a, np.radians(pd))) * af_f
        mb = (1 - sine_mod(l_fine, amp_b, np.radians(pd))) * an_f
        lw = 2.8 if abs(pd - 20) < 0.5 else 1.4
        lab = rf"$\phi={pd}^\circ$" + (" (sheet)" if abs(pd - 20) < 0.5 else "")
        col = "#E69F00" if abs(pd - 20) < 0.5 else c
        ax0.plot(l_fine, ma, color=col, lw=lw, label=lab, zorder=3 if abs(pd-20)<0.5 else 2)
        ax1.plot(l_fine, mb, color=col, lw=lw, zorder=3 if abs(pd-20)<0.5 else 2)

    ax0.legend(loc="upper right", fontsize=7.5, ncol=2, framealpha=0.92)
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.92)
    ax1.set_xlabel(r"Galactic longitude $\ell$ [deg]")
    fig.suptitle(
        r"Migkas Fig.6-style — $\phi$ level sets (peaks slide in $\ell$)" "\n"
        + rf"$\phi$ levels at Amp $=\pm 0.06$ — sheet $\phi\approx 20^\circ$: "
          rf"$R^2_A={r2_a_sheet:.3f}$, $R^2_B={r2_b_sheet:.3f}$",
        fontsize=12, y=0.98,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"levelsets_phi.{ext}", dpi=300)
    plt.close(fig)
    print("Wrote levelsets_phi")

    # --- Amp level sets at φ ≈ 20° ---
    amps = [0.00, -0.03, -0.06, -0.09, -0.12]
    colors = [cmap(x) for x in np.linspace(0.1, 0.92, len(amps))]

    fig = plt.figure(figsize=(8.2, 7.2))
    gs = GridSpec(2, 1, height_ratios=[1, 1], hspace=0.08)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    for ax, y, elo, ehi, ylab, which in (
        (ax0, mig["A"], mig["A_elo"], mig["A_ehi"], r"$A/A_{\mathrm{all}}$", "A"),
        (ax1, mig["B"], mig["B_elo"], mig["B_ehi"], r"$B/B_{\mathrm{all}}$", "B"),
    ):
        ax.errorbar(
            l, y, yerr=[elo, ehi], fmt="o", ms=3.5, color="k",
            ecolor="#999999", elinewidth=0.8, capsize=0, zorder=5,
            label="Migkas" if which == "A" else r"Migkas $B$",
        )
        ax.axhline(1.0, color="#888888", ls="--", lw=0.9)
        ax.set_ylabel(ylab)
        ax.set_ylim(0.80, 1.30)
        ax.grid(True, alpha=0.35)
        ax.set_xlim(0, 360)

    for amp, c in zip(amps, colors):
        ma = (1 - sine_mod(l_fine, amp, SHEET_PHI_RAD)) * af_f
        mb = (1 - sine_mod(l_fine, -amp, SHEET_PHI_RAD)) * an_f
        is_sheet = abs(amp - SHEET_AMP_A) < 1e-9
        lw = 2.8 if is_sheet else 1.4
        lab = rf"Amp$_A={amp:+.2f}$" + (" (sheet)" if is_sheet else "")
        col = "#E69F00" if is_sheet else c
        ax0.plot(l_fine, ma, color=col, lw=lw, label=lab)
        ax1.plot(l_fine, mb, color=col, lw=lw)

    ax0.legend(loc="upper right", fontsize=8, framealpha=0.92)
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.92)
    ax1.set_xlabel(r"Galactic longitude $\ell$ [deg]")
    fig.suptitle(
        r"Migkas Fig.6-style — Amp level sets (features grow / shrink)" "\n"
        + rf"Amp levels at $\phi\approx 20^\circ$ — modulation depth",
        fontsize=12, y=0.98,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"levelsets_Amp.{ext}", dpi=300)
    plt.close(fig)
    print("Wrote levelsets_Amp")

    return dict(R2_A_sheet=r2_a_sheet, R2_B_sheet=r2_b_sheet)


# -------------------- Part 2: geometry --------------------

def scale_kappa(pmeta: dict, kappa_target: float) -> dict:
    """Scale masses so κ = GM_enc/(D³H²) hits target (breaks ρ(D)=Ω_mρ_c)."""
    k0 = float(pmeta["kappa"])
    if k0 <= 0:
        raise ValueError("kappa<=0")
    s = float(kappa_target) / k0
    p = dict(pmeta)
    for key in ("M_tot", "M_enc", "kappa", "kappa_core", "kappa_u", "rho0", "rho_u",
                "rho_local_over_rcrit", "rho_mean_over_rcrit"):
        if key in p and p[key] is not None:
            p[key] = float(p[key]) * s
    return p


def af_an_alpha(pmeta, eps, alpha, af_mode="hratio", n_phi=361):
    """AF/AN with tunable α (velocity scale) and eps=d/D."""
    l_deg = np.linspace(0.0, 360.0, n_phi, endpoint=False)
    phi = np.deg2rad(l_deg)
    n_hat = np.stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)], axis=-1)
    d = eps * D
    D_vec = np.array([D, 0.0, 0.0])
    r_vec = D_vec + d * n_hat
    u_r = peculiar_u(r_vec, pmeta, alpha=alpha)
    u_D = peculiar_u(D_vec.reshape(1, 3), pmeta, alpha=alpha)[0]
    delta_v = np.sum((u_r - u_D) * n_hat, axis=-1)
    H_obs = H + delta_v / d
    m_src = m_enc_at_r(np.linalg.norm(r_vec, axis=-1), pmeta)
    H_safe = np.maximum(H_obs, 1e-6 * H)
    if af_mode == "invsq":
        A_raw = 1.0 / (H_safe / H) ** 2
    else:
        A_raw = H_safe / H
    AF = A_raw / np.mean(A_raw)
    AN = (m_src / pmeta["M_enc"]) / np.mean(m_src / pmeta["M_enc"])
    return l_deg, boxcar_dl90(AF, l_deg), boxcar_dl90(AN, l_deg)


def model_AB_fixed(af, an, l, amp_a=SHEET_AMP_A, amp_b=BEST_AMP_B,
                   phi_a_deg=SHEET_PHI_DEG, phi_b_deg=BEST_PHI_B_DEG):
    ma = (1.0 - sine_mod(l, amp_a, np.radians(phi_a_deg))) * af
    mb = (1.0 - sine_mod(l, amp_b, np.radians(phi_b_deg))) * an
    return ma, mb


def r2_geometry_cell(mig, x=None, f_u=None, alpha=None, eps=None, kappa=None,
                     af_mode="hratio"):
    """Build pmeta, optionally rescale κ, return R²_A, R²_B with fixed Amp×φ."""
    x = BASE["x"] if x is None else x
    f_u = BASE["f_u"] if f_u is None else f_u
    alpha = BASE["alpha"] if alpha is None else alpha
    eps = BASE["eps"] if eps is None else eps
    pmeta = params_twoscale_local_match(float(x), float(f_u))
    if kappa is not None:
        pmeta = scale_kappa(pmeta, float(kappa))
    l_deg, AF, AN = af_an_alpha(pmeta, eps=float(eps), alpha=float(alpha), af_mode=af_mode)
    af_m = interp_periodic(mig["l"], l_deg, AF)
    an_m = interp_periodic(mig["l"], l_deg, AN)
    ma, mb = model_AB_fixed(af_m, an_m, mig["l"])
    return dict(
        R2_A=rsq(mig["A"], ma),
        R2_B=rsq(mig["B"], mb),
        kappa=float(pmeta["kappa"]),
        af_m=af_m,
        an_m=an_m,
        l_deg=l_deg,
        AF=AF,
        AN=AN,
    )


def heat_ab(xv, yv, ZA, ZB, xlab, ylab, title, name, marks=(), hlines=(), vlines=(),
            out_prefix="heatmap"):
    """Write separate A and B magma/plasma heatmaps (Amp×φ style)."""
    meta = {}
    for Z, tag, cmap, ylab_cb in (
        (ZA, "A", "magma", r"$R^2_A$ vs Migkas $A/A_{\mathrm{all}}$"),
        (ZB, "B", "plasma", r"$R^2_B$ vs Migkas $B/B_{\mathrm{all}}$"),
    ):
        fig, ax = plt.subplots(figsize=(6.8, 4.8))
        pm = ax.pcolormesh(
            edges(xv), edges(yv), Z, cmap=cmap, shading="flat", vmin=0.0,
            edgecolors="face", linewidth=0, rasterized=True,
        )
        iy, ix = np.unravel_index(np.nanargmax(Z), Z.shape)
        for y0, lab in hlines:
            ax.axhline(y0, color="#17becf", ls=":", lw=1.5, label=lab)
        for x0, lab in vlines:
            ax.axvline(x0, color="#17becf", ls=":", lw=1.5, label=lab)
        for (mx, my, mk, lab) in marks:
            ax.plot(mx, my, marker=mk, ms=10, color="#17becf", mec="k", mew=0.7,
                    ls="none", label=lab)
        ax.plot(
            xv[ix], yv[iy], marker="*", ms=16, color="w", mec="k", mew=0.8,
            ls="none",
            label=rf"best: ({xv[ix]:.3g}, {yv[iy]:.3g}), $R^2={Z[iy, ix]:.3f}$",
        )
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(f"{title} — $R^2_{tag}$", fontsize=11)
        cb = fig.colorbar(pm, ax=ax)
        cb.set_label(ylab_cb)
        ax.legend(loc="best", fontsize=7.5, framealpha=0.9)
        fig.tight_layout()
        stem = f"{out_prefix}_{name}_{tag}"
        for ext in ("png", "pdf"):
            fig.savefig(OUT / f"{stem}.{ext}", dpi=300)
        plt.close(fig)
        meta[tag] = dict(best_x=float(xv[ix]), best_y=float(yv[iy]),
                         R2=float(Z[iy, ix]), path=str(OUT / f"{stem}.png"))
        print(f"Wrote {stem}.png  best R2={Z[iy,ix]:.4f}")
    return meta


def levelset_geometry_curves(mig, vary_name, vary_vals, hold: dict, label_fmt,
                             out_stem, title, subtitle):
    """A/B vs ℓ level-set curves while slicing one geometry param."""
    l = mig["l"]
    l_fine = np.linspace(0, 360, 361, endpoint=False)
    cmap = plt.cm.viridis
    colors = [cmap(x) for x in np.linspace(0.08, 0.92, len(vary_vals))]

    fig = plt.figure(figsize=(8.2, 7.2))
    gs = GridSpec(2, 1, height_ratios=[1, 1], hspace=0.08)
    ax0 = fig.add_subplot(gs[0])
    ax1 = fig.add_subplot(gs[1], sharex=ax0)

    for ax, y, elo, ehi, ylab, which in (
        (ax0, mig["A"], mig["A_elo"], mig["A_ehi"], r"$A/A_{\mathrm{all}}$", "A"),
        (ax1, mig["B"], mig["B_elo"], mig["B_ehi"], r"$B/B_{\mathrm{all}}$", "B"),
    ):
        ax.errorbar(
            l, y, yerr=[elo, ehi], fmt="o", ms=3.5, color="k",
            ecolor="#999999", elinewidth=0.8, capsize=0, zorder=5,
            label="Migkas" if which == "A" else r"Migkas $B/B_{\mathrm{all}}$",
        )
        ax.axhline(1.0, color="#888888", ls="--", lw=0.9)
        ax.set_ylabel(ylab)
        ax.set_ylim(0.75, 1.35)
        ax.grid(True, alpha=0.35)
        ax.set_xlim(0, 360)

    sheet_key = None
    for val, c in zip(vary_vals, colors):
        kw = dict(hold)
        kw[vary_name] = val
        # resolve baseline defaults
        rec = r2_geometry_cell(
            mig,
            x=kw.get("x", BASE["x"]),
            f_u=kw.get("f_u", BASE["f_u"]),
            alpha=kw.get("alpha", BASE["alpha"]),
            eps=kw.get("eps", BASE["eps"]),
            kappa=kw.get("kappa", None),
        )
        # rebuild fine curves
        pmeta = params_twoscale_local_match(float(kw.get("x", BASE["x"])),
                                            float(kw.get("f_u", BASE["f_u"])))
        if kw.get("kappa") is not None:
            pmeta = scale_kappa(pmeta, float(kw["kappa"]))
        lf, AF, AN = af_an_alpha(
            pmeta, eps=float(kw.get("eps", BASE["eps"])),
            alpha=float(kw.get("alpha", BASE["alpha"])),
        )
        ma, mb = model_AB_fixed(AF, AN, lf)
        is_base = False
        if vary_name == "eps" and abs(val - BASE["eps"]) < 1e-9:
            is_base = True
        if vary_name == "alpha" and abs(val - BASE["alpha"]) < 1e-9:
            is_base = True
        if vary_name == "x" and abs(val - BASE["x"]) < 1e-9:
            is_base = True
        if vary_name == "kappa" and abs(val - 0.2306) < 0.01:
            is_base = True
        lw = 2.6 if is_base else 1.4
        col = "#E69F00" if is_base else c
        lab = label_fmt(val) + (" (sheet/best)" if is_base else "")
        ax0.plot(lf, ma, color=col, lw=lw, label=lab)
        ax1.plot(lf, mb, color=col, lw=lw)

    ax0.legend(loc="upper right", fontsize=7.5, ncol=2, framealpha=0.92)
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.92)
    ax1.set_xlabel(r"Galactic longitude $\ell$ [deg]")
    fig.suptitle(title + "\n" + subtitle, fontsize=11.5, y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{out_stem}.{ext}", dpi=300)
    plt.close(fig)
    print(f"Wrote {out_stem}")


def part2_geometry(mig):
    t0 = time.time()
    report = {}

    # Pair 1: κ (mass) × d/D (obs depth / R)
    kappa_g = np.round(np.arange(0.05, 1.0001, 0.05), 4)
    eps_g = np.round(np.arange(0.02, 0.5001, 0.02), 3)
    ZA = np.zeros((len(eps_g), len(kappa_g)))
    ZB = np.zeros_like(ZA)
    print(f"Scanning κ × d/D  ({len(kappa_g)}×{len(eps_g)}) ...", flush=True)
    for j, eps in enumerate(eps_g):
        for i, kap in enumerate(kappa_g):
            rec = r2_geometry_cell(mig, kappa=kap, eps=eps)
            ZA[j, i], ZB[j, i] = rec["R2_A"], rec["R2_B"]
        if (j + 1) % 5 == 0:
            print(f"  eps row {j+1}/{len(eps_g)}  ({time.time()-t0:.1f}s)", flush=True)
    report["kappa_vs_dD"] = heat_ab(
        kappa_g, eps_g, ZA, ZB,
        r"$\kappa = GM(<D)/(D^3 H^2)$ (mass)",
        r"$d/D$ (obs. depth / vision radius)",
        r"Geometry pair: mass $\kappa$ × depth $d/D$" "\n"
        rf"$D/L={BASE['x']}$, $f_u={BASE['f_u']}$, $\alpha={BASE['alpha']}$, "
        r"AF=hratio, fixed Amp×φ (sheet A / best B)",
        "kappa_dD",
        marks=[(0.2306, 0.10, "o", r"baseline $\kappa\approx0.231$, $d/D=0.1$")],
        vlines=[(0.15, r"$\kappa=\Omega_m/2=0.15$")],
    )
    report["kappa_vs_dD"]["ZA_max"] = float(np.nanmax(ZA))
    report["kappa_vs_dD"]["ZB_max"] = float(np.nanmax(ZB))

    # Pair 2: D/L (distance scale) × α (velocity scale)
    x_g = np.round(np.arange(3.0, 10.001, 0.25), 3)
    alpha_g = np.round(np.arange(0.2, 2.001, 0.1), 3)
    ZA2 = np.zeros((len(alpha_g), len(x_g)))
    ZB2 = np.zeros_like(ZA2)
    print(f"Scanning D/L × α  ({len(x_g)}×{len(alpha_g)}) ...", flush=True)
    for j, a in enumerate(alpha_g):
        for i, x in enumerate(x_g):
            rec = r2_geometry_cell(mig, x=x, alpha=a)
            ZA2[j, i], ZB2[j, i] = rec["R2_A"], rec["R2_B"]
        if (j + 1) % 4 == 0:
            print(f"  alpha row {j+1}/{len(alpha_g)}  ({time.time()-t0:.1f}s)", flush=True)
    report["DL_vs_alpha"] = heat_ab(
        x_g, alpha_g, ZA2, ZB2,
        r"$D/L$ (CoE→observer distance / scale length)",
        r"$\alpha$ (velocity scale in $u=-\alpha\,GM/(r^2 H)\,\hat r$)",
        r"Geometry pair: $D/L$ × velocity $\alpha$" "\n"
        rf"$f_u={BASE['f_u']}$, $d/D={BASE['eps']}$, AF=hratio, "
        r"κ from local-match, fixed Amp×φ",
        "DL_alpha",
        marks=[(BASE["x"], BASE["alpha"], "o", r"baseline $D/L=5.9$, $\alpha=1$")],
    )
    report["DL_vs_alpha"]["ZA_max"] = float(np.nanmax(ZA2))
    report["DL_vs_alpha"]["ZB_max"] = float(np.nanmax(ZB2))

    # Pair 3 (bonus): D/L × d/D — distance × depth
    ZA3 = np.zeros((len(eps_g), len(x_g)))
    ZB3 = np.zeros_like(ZA3)
    print(f"Scanning D/L × d/D  ({len(x_g)}×{len(eps_g)}) ...", flush=True)
    for j, eps in enumerate(eps_g):
        for i, x in enumerate(x_g):
            rec = r2_geometry_cell(mig, x=x, eps=eps)
            ZA3[j, i], ZB3[j, i] = rec["R2_A"], rec["R2_B"]
    report["DL_vs_dD"] = heat_ab(
        x_g, eps_g, ZA3, ZB3,
        r"$D/L$",
        r"$d/D$",
        r"Geometry pair: $D/L$ × $d/D$" "\n"
        rf"$\alpha={BASE['alpha']}$, $f_u={BASE['f_u']}$, local-match κ, fixed Amp×φ",
        "DL_dD",
        marks=[(BASE["x"], BASE["eps"], "o", "baseline")],
    )

    # Level-set curves for each dial
    levelset_geometry_curves(
        mig, "kappa",
        [0.10, 0.15, 0.23, 0.40, 0.60, 0.90],
        dict(eps=BASE["eps"], alpha=BASE["alpha"], x=BASE["x"]),
        lambda v: rf"$\kappa={v:.2f}$",
        "levelset_kappa",
        r"Migkas Fig.6-style — mass $\kappa$ level sets (features grow / shrink)",
        rf"$\kappa$ slices at $d/D={BASE['eps']}$, $D/L={BASE['x']}$, $\alpha={BASE['alpha']}$",
    )
    levelset_geometry_curves(
        mig, "eps",
        [0.02, 0.05, 0.10, 0.15, 0.25, 0.40],
        dict(kappa=0.2306, alpha=BASE["alpha"], x=BASE["x"]),
        lambda v: rf"$d/D={v:.2f}$",
        "levelset_dD",
        r"Migkas Fig.6-style — $d/D$ (obs. depth) level sets",
        rf"$d/D$ slices at $\kappa\approx0.231$, $D/L={BASE['x']}$, $\alpha={BASE['alpha']}$",
    )
    levelset_geometry_curves(
        mig, "alpha",
        [0.3, 0.6, 1.0, 1.4, 1.8],
        dict(eps=BASE["eps"], x=BASE["x"]),  # κ from local-match
        lambda v: rf"$\alpha={v:.1f}$",
        "levelset_alpha",
        r"Migkas Fig.6-style — velocity scale $\alpha$ level sets",
        rf"$\alpha$ slices at $d/D={BASE['eps']}$, $D/L={BASE['x']}$ (κ from local-match)",
    )
    levelset_geometry_curves(
        mig, "x",
        [3.5, 4.5, 5.9, 7.5, 9.0],
        dict(eps=BASE["eps"], alpha=BASE["alpha"]),
        lambda v: rf"$D/L={v:.1f}$",
        "levelset_DL",
        r"Migkas Fig.6-style — $D/L$ (distance scale) level sets",
        rf"$D/L$ slices at $d/D={BASE['eps']}$, $\alpha={BASE['alpha']}$ (κ from local-match)",
    )

    report["runtime_s"] = time.time() - t0
    return report


def main():
    style()
    OUT.mkdir(parents=True, exist_ok=True)
    mig = load_migkas()
    sheet = load_sheet()
    print(f"Migkas n={len(mig['l'])}, sheet n={len(sheet['AR'])}")

    print("\n=== Part 1: Amp × φ (sheet AF/AN) ===")
    p1 = part1_heatmaps(mig, sheet)
    p1_ls = part1_levelsets(mig, sheet, p1["af"], p1["an"])

    print("\n=== Part 2: geometry pairs ===")
    p2 = part2_geometry(mig)

    meta = {
        "part1": {
            "model": "BA=Amp*sin(2πℓ/360+φ); A=(1-BA)*AF; B=(1-Amp_B*sin)*AN; Pearson R²",
            "AF_AN_source": str(SHEET_CSV),
            "migkas": str(MIG_CSV),
            "best_A": p1["A"],
            "best_B": p1["B"],
            "sheet_point": {
                "Amp": SHEET_AMP_A, "phi_rad": SHEET_PHI_RAD,
                "phi_deg": SHEET_PHI_DEG,
                "R2_A": p1["sheet_R2_A"],
                "R2_B_opp": p1_ls["R2_B_sheet"],
            },
            "outputs": [
                str(OUT / "heatmap_A_Amp_phi.png"),
                str(OUT / "heatmap_B_AmpB_phi.png"),
                str(OUT / "levelsets_phi.png"),
                str(OUT / "levelsets_Amp.png"),
            ],
        },
        "part2": {
            "baseline": BASE,
            "pairs": ["kappa_vs_dD", "DL_vs_alpha", "DL_vs_dD"],
            "covers": {
                "D": "D/L in DL_vs_alpha, DL_vs_dD, levelset_DL",
                "R_obs_depth": "d/D in kappa_vs_dD, DL_vs_dD, levelset_dD",
                "velocity": "alpha in DL_vs_alpha, levelset_alpha",
                "mass": "kappa in kappa_vs_dD, levelset_kappa",
            },
            "results": {k: {kk: vv for kk, vv in v.items() if kk not in ("ZA", "ZB")}
                        if isinstance(v, dict) else v
                        for k, v in p2.items()},
            "outputs_glob": str(OUT / "heatmap_*_{A,B}.png") + ", " + str(OUT / "levelset_*.png"),
        },
    }
    # strip non-jsonable
    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()
                    if not isinstance(v, np.ndarray)}
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        return o

    meta = clean(meta)
    json.dump(meta, open(OUT / "heatmap_geometry_report.json", "w"), indent=2)
    print("\n=== DONE ===")
    print(json.dumps(meta, indent=2)[:2500])
    return meta


if __name__ == "__main__":
    main()

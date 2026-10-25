#!/usr/bin/env python3
"""
Level-set (contour) maps of R²_A and R²_B vs model parameters.

Baseline: two-scale ST+Gauss + Migkas under ρ(D)=Ω_m ρ_crit=0.3 ρ_crit
  D/L≈5.9, f_u≈0.94, κ≈0.23, AF=hratio, Amp=+0.07, φ≈213°, α=1, d/D=0.1

Reuses physics from run_density_local_match.py and Amp/φ fitters from
r2_maximize/run_r2_maximize.py (independent Amp_A/Amp_B, shared φ).
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent  # density_local_match/
ROOT = Path("/workspace/cosmology")
OUT = HERE

# ---- load density_local_match helpers ----
_spec_d = importlib.util.spec_from_file_location(
    "run_density_local_match", PARENT / "run_density_local_match.py"
)
_dmod = importlib.util.module_from_spec(_spec_d)
sys.modules["run_density_local_match"] = _dmod
_spec_d.loader.exec_module(_dmod)

# ---- load r2_maximize helpers (alpha/eps-aware physics + indep Amp fit) ----
_spec_r = importlib.util.spec_from_file_location(
    "run_r2_maximize", PARENT / "r2_maximize" / "run_r2_maximize.py"
)
_rmod = importlib.util.module_from_spec(_spec_r)
sys.modules["run_r2_maximize"] = _rmod
_spec_r.loader.exec_module(_rmod)

OMEGA_M = _dmod.OMEGA_M
OMEGA = _dmod.OMEGA
C_BLUE = _dmod.C_BLUE
C_ORANGE = _dmod.C_ORANGE
C_GREEN = _dmod.C_GREEN
C_BLACK = _dmod.C_BLACK
C_RED = _dmod.C_RED
C_PURPLE = _dmod.C_PURPLE
style = _dmod.style
load_migkas = _dmod.load_migkas
interp_periodic = _dmod.interp_periodic
rsq = _dmod.rsq
sine_mod = _dmod.sine_mod

params_twoscale_local_match = _rmod.params_twoscale_local_match
make_pmeta = _rmod.make_pmeta
af_an_from_physics = _rmod.af_an_from_physics
fit_indep_amps_shared_phi = _rmod.fit_indep_amps_shared_phi
fit_shared_opp = _rmod.fit_shared_opp
_amp_grid = _rmod._amp_grid
model_AB = _rmod.model_AB
model_A_only = _rmod.model_A_only
model_B_only = _rmod.model_B_only

# ---- baseline (balanced recommended_profile.json) ----
BASE = dict(
    x=5.9,
    f_u=0.94,
    kappa=0.23058828928336284,
    af_mode="hratio",
    amp=0.07,
    phi=3.7143692147415237,
    phi_deg=212.81767955801106,
    alpha=1.0,
    eps=0.1,
    R2_A=0.6852166343310685,
    R2_B=0.6854567749942476,
)

# Coarse but readable grids; light Amp/φ re-fit per cell
N_PHI_FIT = 73          # ~5° steps
AMP_GRID_FIT = _amp_grid(0.02, 0.22, 21)
N_PHYS_PHI = 361         # sky circle resolution (lighter than 721)


def style_levelsets():
    style()
    mpl.rcParams.update(
        {
            "savefig.dpi": 200,
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 11,
            "legend.fontsize": 8,
        }
    )


def physics_templates(pmeta, alpha, eps, af_mode, mig):
    """AF/AN on Migkas ℓ grid + κ check."""
    assert abs(pmeta["rho_local_over_rcrit"] - OMEGA_M) < 0.01 * OMEGA_M + 1e-9
    phys = af_an_from_physics(pmeta, eps=eps, alpha=alpha, af_mode=af_mode)
    # recompute with lighter n_phi if needed — af_an uses default 721; OK
    af_m = interp_periodic(mig["l"], phys["l_deg"], phys["AF_smooth"])
    an_m = interp_periodic(mig["l"], phys["l_deg"], phys["AN_smooth"])
    return af_m, an_m, float(pmeta["kappa"]), float(phys["frac_ptp"])


def refit_amp_phi(af_m, an_m, mig):
    """Independent Amp_A, Amp_B; shared φ; balanced = max min(R²_A,R²_B)."""
    indep = fit_indep_amps_shared_phi(
        mig["l"], af_m, an_m, mig["A"], mig["B"],
        amp_grid=AMP_GRID_FIT, n_phi=N_PHI_FIT,
    )
    shared = fit_shared_opp(
        mig["l"], af_m, an_m, mig["A"], mig["B"],
        amp_grid=AMP_GRID_FIT, n_phi=N_PHI_FIT,
    )
    cand = [indep["best_min"], shared["best_min"]]
    best = max(cand, key=lambda r: (r["r2_min"], r["r2_sum"]))
    return best


def eval_cell(x, f_u, alpha, eps, af_mode, mig):
    pmeta = make_pmeta("twoscale", x, f_u=f_u)
    af_m, an_m, kappa, frac_ptp = physics_templates(pmeta, alpha, eps, af_mode, mig)
    fit = refit_amp_phi(af_m, an_m, mig)
    return {
        "x": float(x),
        "f_u": float(f_u),
        "alpha": float(alpha),
        "eps": float(eps),
        "kappa": kappa,
        "frac_ptp": frac_ptp,
        "R2_A": fit["r2_a"],
        "R2_B": fit["r2_b"],
        "R2_min": fit["r2_min"],
        "Amp_A": fit["amp_a"],
        "Amp_B": fit["amp_b"],
        "phi": fit["phi"],
        "phi_deg": fit["phi_deg"],
        "mode": fit["mode"],
    }


def scan_grid(xs, ys, cell_fn, label=""):
    """xs, ys 1d arrays; cell_fn(i,j,x,y) -> dict with R2_A, R2_B, kappa."""
    nx, ny = len(xs), len(ys)
    R2A = np.full((ny, nx), np.nan)
    R2B = np.full((ny, nx), np.nan)
    KAP = np.full((ny, nx), np.nan)
    rows = []
    t0 = time.time()
    total = nx * ny
    done = 0
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            rec = cell_fn(i, j, x, y)
            R2A[j, i] = rec["R2_A"]
            R2B[j, i] = rec["R2_B"]
            KAP[j, i] = rec.get("kappa", np.nan)
            rows.append(rec)
            done += 1
            if done % 20 == 0 or done == total:
                dt = time.time() - t0
                print(
                    f"  [{label}] {done}/{total}  "
                    f"({dt:.1f}s, {dt/max(done,1):.2f}s/cell)",
                    flush=True,
                )
    return dict(xs=xs, ys=ys, R2_A=R2A, R2_B=R2B, kappa=KAP, rows=rows)


def _levels(Z, n=9):
    z = Z[np.isfinite(Z)]
    if z.size == 0:
        return np.linspace(0, 1, n)
    lo, hi = float(np.nanpercentile(z, 5)), float(np.nanpercentile(z, 99))
    if hi - lo < 1e-3:
        lo, hi = float(np.nanmin(z)), float(np.nanmax(z)) + 1e-6
    return np.linspace(lo, hi, n)


def plot_ab_pair(
    xs, ys, R2A, R2B, kappa,
    xlabel, ylabel, title,
    baseline_xy, out_stem,
    kappa_overlay=True,
    xscale=None, yscale=None,
):
    """Side-by-side R²_A / R²_B filled contours; optional κ dashed overlay."""
    style_levelsets()
    X, Y = np.meshgrid(xs, ys)
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), constrained_layout=True)

    for ax, Z, name, cmap in (
        (axes[0], R2A, r"$R^2_A$", "viridis"),
        (axes[1], R2B, r"$R^2_B$", "plasma"),
    ):
        lev = _levels(Z, n=10)
        cf = ax.contourf(X, Y, Z, levels=lev, cmap=cmap, extend="both")
        cs = ax.contour(X, Y, Z, levels=lev, colors="k", linewidths=0.35, alpha=0.45)
        ax.clabel(cs, fmt="%.2f", fontsize=7, inline=True)
        cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(name)

        if kappa_overlay and np.isfinite(kappa).any():
            klev = np.unique(np.round(_levels(kappa, n=6), 3))
            if klev.size >= 2:
                ck = ax.contour(
                    X, Y, kappa, levels=klev, colors="white",
                    linewidths=1.1, linestyles="--", alpha=0.85,
                )
                ax.clabel(ck, fmt=r"$\kappa$=%.2f", fontsize=7)

        if baseline_xy is not None:
            ax.plot(
                baseline_xy[0], baseline_xy[1],
                marker="*", markersize=16, color="white",
                markeredgecolor="black", markeredgewidth=0.9,
                zorder=10, label="baseline",
            )
            ax.legend(loc="upper right", framealpha=0.9, fontsize=8)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{name}  —  {title}")
        if xscale:
            ax.set_xscale(xscale)
        if yscale:
            ax.set_yscale(yscale)

    fig.suptitle(
        r"ST+Gauss + Migkas  |  $\rho(D)=\Omega_m\rho_{\rm crit}$  |  AF=hratio"
        + (r"  |  Amp$_A$, Amp$_B$, $\phi$ re-fit per cell" if "Amp" not in title else ""),
        fontsize=10, y=1.02,
    )
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{out_stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_stem}.png/pdf", flush=True)


def plot_amp_phi_landscape(af_m, an_m, mig, out_stem="levelsets_Amp_vs_phi"):
    """Diagnostic Amp×φ landscape for A and for B separately at fixed baseline physics."""
    style_levelsets()
    amps = np.linspace(-0.25, 0.25, 41)
    phis_deg = np.linspace(0.0, 360.0, 73, endpoint=False)
    phis = np.deg2rad(phis_deg)
    R2A = np.zeros((len(amps), len(phis)))
    R2B = np.zeros((len(amps), len(phis)))
    l = mig["l"]
    for i, amp in enumerate(amps):
        for j, phi in enumerate(phis):
            ma = model_A_only(af_m, l, amp, phi)
            mb = model_B_only(an_m, l, amp, phi)
            R2A[i, j] = rsq(mig["A"], ma)
            R2B[i, j] = rsq(mig["B"], mb)

    PHI, AMP = np.meshgrid(phis_deg, amps)
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6), constrained_layout=True)
    for ax, Z, name, cmap in (
        (axes[0], R2A, r"$R^2_A$(Amp$_A$, $\phi$)", "viridis"),
        (axes[1], R2B, r"$R^2_B$(Amp$_B$, $\phi$)", "plasma"),
    ):
        lev = _levels(Z, n=12)
        cf = ax.contourf(PHI, AMP, Z, levels=lev, cmap=cmap, extend="both")
        cs = ax.contour(PHI, AMP, Z, levels=lev, colors="k", linewidths=0.35, alpha=0.4)
        ax.clabel(cs, fmt="%.2f", fontsize=7)
        cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(r"$R^2$")
        # baseline Amp=+0.07, φ≈213° (shared-opp balanced uses Amp_B=-Amp_A;
        # for solo landscape mark Amp=+0.07 at baseline φ)
        ax.plot(
            BASE["phi_deg"], BASE["amp"],
            marker="*", markersize=16, color="white",
            markeredgecolor="black", markeredgewidth=0.9,
            zorder=10, label="baseline Amp,+φ",
        )
        # also mark Amp=-0.07 (opposite) for shared-opp B convention
        ax.plot(
            BASE["phi_deg"], -BASE["amp"],
            marker="o", markersize=7, color="white",
            markeredgecolor="black", markeredgewidth=0.8,
            zorder=9, label="−Amp @ same φ",
        )
        ax.axhline(0, color="0.5", lw=0.6, alpha=0.6)
        ax.set_xlabel(r"$\phi$ [deg]")
        ax.set_ylabel("Amp")
        ax.set_title(name)
        ax.legend(loc="lower right", fontsize=7, framealpha=0.9)
        ax.set_xlim(0, 360)

    fig.suptitle(
        r"Amp–$\phi$ landscape at fixed baseline physics "
        r"(D/L=5.9, $f_u$=0.94, $\alpha$=1, $d/D$=0.1, AF=hratio)",
        fontsize=10, y=1.02,
    )
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{out_stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_stem}.png/pdf", flush=True)
    return dict(amps=amps, phis_deg=phis_deg, R2_A=R2A, R2_B=R2B)


def plot_summary(grids, amp_phi, out_stem="levelsets_summary"):
    """Combined multi-panel figure + annotation strip."""
    style_levelsets()
    fig = plt.figure(figsize=(13.5, 11.5))
    gs = GridSpec(
        3, 4, figure=fig,
        height_ratios=[1.0, 1.0, 0.42],
        hspace=0.38, wspace=0.32,
        left=0.07, right=0.98, top=0.93, bottom=0.04,
    )

    panels = [
        ("DL_fu", r"$D/L$", r"$f_u$", "D/L vs $f_u$", True),
        ("alpha_eps", r"$\alpha$", r"$d/D$", r"$\alpha$ vs $d/D$", False),
        ("DL_alpha", r"$D/L$", r"$\alpha$", r"$D/L$ vs $\alpha$", True),
    ]

    for col, (key, xlab, ylab, ttl, kover) in enumerate(panels):
        g = grids[key]
        X, Y = np.meshgrid(g["xs"], g["ys"])
        for row, (Z, name, cmap) in enumerate(
            ((g["R2_A"], r"$R^2_A$", "viridis"), (g["R2_B"], r"$R^2_B$", "plasma"))
        ):
            ax = fig.add_subplot(gs[row, col])
            lev = _levels(Z, n=8)
            cf = ax.contourf(X, Y, Z, levels=lev, cmap=cmap, extend="both")
            ax.contour(X, Y, Z, levels=lev, colors="k", linewidths=0.25, alpha=0.35)
            if kover and np.isfinite(g["kappa"]).any():
                klev = np.unique(np.round(_levels(g["kappa"], n=5), 3))
                if klev.size >= 2:
                    ax.contour(
                        X, Y, g["kappa"], levels=klev, colors="white",
                        linewidths=0.9, linestyles="--", alpha=0.8,
                    )
            bx, by = g["baseline_xy"]
            ax.plot(bx, by, "*", ms=12, color="white",
                    markeredgecolor="k", markeredgewidth=0.7, zorder=10)
            ax.set_xlabel(xlab, fontsize=9)
            ax.set_ylabel(ylab, fontsize=9)
            ax.set_title(f"{name}: {ttl}", fontsize=9)
            fig.colorbar(cf, ax=ax, fraction=0.05, pad=0.02)

    # Amp–φ landscapes in column 4
    PHI, AMP = np.meshgrid(amp_phi["phis_deg"], amp_phi["amps"])
    for row, (Z, name, cmap) in enumerate(
        ((amp_phi["R2_A"], r"$R^2_A$", "viridis"),
         (amp_phi["R2_B"], r"$R^2_B$", "plasma"))
    ):
        ax = fig.add_subplot(gs[row, 3])
        lev = _levels(Z, n=8)
        cf = ax.contourf(PHI, AMP, Z, levels=lev, cmap=cmap, extend="both")
        ax.plot(BASE["phi_deg"], BASE["amp"], "*", ms=12, color="white",
                markeredgecolor="k", markeredgewidth=0.7, zorder=10)
        ax.set_xlabel(r"$\phi$ [deg]", fontsize=9)
        ax.set_ylabel("Amp", fontsize=9)
        ax.set_title(f"{name}: Amp vs φ", fontsize=9)
        ax.set_xlim(0, 360)
        fig.colorbar(cf, ax=ax, fraction=0.05, pad=0.02)

    # annotation strip
    axn = fig.add_subplot(gs[2, :])
    axn.axis("off")
    note = (
        r"$\mathbf{What\ each\ dial\ does}$  (baseline $\star$: "
        r"$D/L{=}5.9$, $f_u{=}0.94$, $\kappa{\approx}0.23$, Amp${=}{+}0.07$, "
        r"$\phi{\approx}213^\circ$, $\alpha{=}1$, $d/D{=}0.1$; "
        r"$\rho(D){=}0.3\,\rho_{\rm crit}$ always)."
        "\n\n"
        r"$\bullet$ $D/L$ & $f_u$: set the two-scale density shape $\Rightarrow$ "
        r"tidal amplitude (frac ptp) and enclosed-mass $\kappa=GM_{\rm enc}/(D^3 H^2)$ "
        r"(white dashed $=$ $\kappa$ contours)."
        "\n"
        r"$\bullet$ $\alpha$: overall gravity / peculiar-velocity strength "
        r"(scales $H_{\rm obs}$ anisotropy)."
        "\n"
        r"$\bullet$ $d/D$: observer offset $\Rightarrow$ multipole mix in AF/AN "
        r"(dipole vs higher)."
        "\n"
        r"$\bullet$ $\phi$: aims the Amp$\times\sin$ template on the sky (phase / axis).  "
        r"Amp$_A$/Amp$_B$: template depth on A vs B "
        r"(shared $\phi$; re-fit per cell on profile panels)."
    )
    axn.text(
        0.01, 0.95, note, transform=axn.transAxes, va="top", ha="left",
        fontsize=9, family="DejaVu Sans",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="#f7f7f7", edgecolor="#cccccc"),
    )

    fig.suptitle(
        r"Level sets of Migkas $R^2_A$ & $R^2_B$  |  ST+Gauss, AF=hratio, "
        r"$\rho(D)=\Omega_m\rho_{\rm crit}$",
        fontsize=12, y=0.98,
    )
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{out_stem}.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_stem}.png/pdf", flush=True)


def write_note(grids, amp_phi, verify):
    path = OUT / "levelsets_NOTE.md"
    g1 = grids["DL_fu"]
    g2 = grids["alpha_eps"]
    g3 = grids["DL_alpha"]
    # peak cells
    def peak(g, key="R2_min"):
        rows = g["rows"]
        if key == "R2_min":
            best = max(rows, key=lambda r: (r["R2_min"], r["R2_A"] + r["R2_B"]))
        else:
            best = max(rows, key=lambda r: r[key])
        return best

    p1 = peak(g1)
    p2 = peak(g2)
    p3 = peak(g3)
    text = f"""# Level-set guide: which dial moves what

Baseline (balanced recommended profile, held when not scanning that axis):

- two-scale, **D/L = 5.9**, **f_u = 0.94**, **κ ≈ 0.231**, AF = `hratio`
- Amp = **+0.07**, φ ≈ **213°** (3.714 rad), **α = 1**, **d/D = 0.1**
- Constraint: **ρ(D) = Ω_m ρ_crit = 0.3 ρ_crit** on every profile cell
- Amp_A, Amp_B independent with **shared φ**, re-fit per cell on profile panels
  (best of indep-amps vs shared-opposite under max-min R²)

## One-sentence dial map

| Dial | What it primarily changes |
|------|---------------------------|
| **D/L** | Core compactness of the exponential core → tidal amplitude of H_obs and **κ** |
| **f_u** | Uniform outer mass fraction → boosts κ with weaker core tide; trades shape |
| **α** | Gravity / peculiar-velocity **strength** (overall anisotropy scale) |
| **d/D (eps)** | Observer offset → **multipole mix** (dipole vs quadrupole in AF/AN) |
| **φ** | Sky **axis aim** of the Amp×sin template |
| **Amp_A / Amp_B** | Template **depth** on A vs B (phase shared) |

κ is reported as white dashed contours wherever the density profile changes
(D/L–f_u and D/L–α). κ is fixed on α–eps and Amp–φ panels (same baseline profile).

## Panel peaks (this run)

- **D/L vs f_u**: best min(R²) ≈ {p1['R2_min']:.3f} at D/L={p1['x']:.2f}, f_u={p1['f_u']:.3f}, κ={p1['kappa']:.3f}
  (R²_A={p1['R2_A']:.3f}, R²_B={p1['R2_B']:.3f}, Amp_A={p1['Amp_A']:+.3f}, Amp_B={p1['Amp_B']:+.3f}, φ={p1['phi_deg']:.1f}°)
- **α vs d/D**: best min(R²) ≈ {p2['R2_min']:.3f} at α={p2['alpha']:.2f}, d/D={p2['eps']:.3f}
  (R²_A={p2['R2_A']:.3f}, R²_B={p2['R2_B']:.3f})
- **D/L vs α**: best min(R²) ≈ {p3['R2_min']:.3f} at D/L={p3['x']:.2f}, α={p3['alpha']:.2f}, κ={p3['kappa']:.3f}

Baseline verification (ρ(D)/ρ_crit): {verify['rho_D']:.4f} (target 0.3).
Baseline re-fit R²_A={verify['R2_A']:.4f}, R²_B={verify['R2_B']:.4f}
(Amp_A={verify['Amp_A']:+.3f}, Amp_B={verify['Amp_B']:+.3f}, φ={verify['phi_deg']:.1f}°).

## Files

- `levelsets_DL_vs_fu.png` (+ pdf)
- `levelsets_alpha_vs_eps.png` (+ pdf)
- `levelsets_Amp_vs_phi.png` (+ pdf)
- `levelsets_DL_vs_alpha.png` (+ pdf)
- `levelsets_summary.png` (+ pdf) — combined + caption strip
- `run_levelsets.py` — regenerable

Digitised Migkas: `/workspace/cosmology/migkas/fig6_digitized.csv` (bars not drawn on contours).
"""
    path.write_text(text)
    print(f"wrote {path}", flush=True)


def main():
    style_levelsets()
    OUT.mkdir(parents=True, exist_ok=True)
    mig = load_migkas()
    af_mode = BASE["af_mode"]

    # ---- verify baseline ----
    p0 = make_pmeta("twoscale", BASE["x"], f_u=BASE["f_u"])
    af0, an0, k0, _ = physics_templates(p0, BASE["alpha"], BASE["eps"], af_mode, mig)
    fit0 = refit_amp_phi(af0, an0, mig)
    verify = {
        "rho_D": float(p0["rho_local_over_rcrit"]),
        "kappa": k0,
        "R2_A": fit0["r2_a"],
        "R2_B": fit0["r2_b"],
        "Amp_A": fit0["amp_a"],
        "Amp_B": fit0["amp_b"],
        "phi_deg": fit0["phi_deg"],
    }
    print(
        f"Baseline verify: ρ(D)/ρc={verify['rho_D']:.4f}, κ={verify['kappa']:.4f}, "
        f"R²_A={verify['R2_A']:.4f}, R²_B={verify['R2_B']:.4f}, "
        f"Amp_A={verify['Amp_A']:+.3f}, Amp_B={verify['Amp_B']:+.3f}, "
        f"φ={verify['phi_deg']:.1f}°",
        flush=True,
    )

    grids = {}

    # ---- 1. D/L vs f_u ----
    xs = np.linspace(3.0, 8.0, 14)       # D/L
    ys = np.linspace(0.85, 0.98, 12)     # f_u
    print("\n=== Panel 1: D/L vs f_u ===", flush=True)

    def cell_dl_fu(i, j, x, y):
        return eval_cell(x, y, BASE["alpha"], BASE["eps"], af_mode, mig)

    g = scan_grid(xs, ys, cell_dl_fu, label="DL_fu")
    g["baseline_xy"] = (BASE["x"], BASE["f_u"])
    grids["DL_fu"] = g
    plot_ab_pair(
        xs, ys, g["R2_A"], g["R2_B"], g["kappa"],
        xlabel=r"$D/L$", ylabel=r"$f_u$",
        title=r"$D/L$ vs $f_u$ (α=1, $d/D$=0.1)",
        baseline_xy=g["baseline_xy"],
        out_stem="levelsets_DL_vs_fu",
        kappa_overlay=True,
    )

    # ---- 2. α vs eps ----
    xs2 = np.linspace(0.5, 2.0, 14)
    ys2 = np.linspace(0.04, 0.22, 12)
    print("\n=== Panel 2: α vs d/D ===", flush=True)

    def cell_a_e(i, j, a, e):
        return eval_cell(BASE["x"], BASE["f_u"], a, e, af_mode, mig)

    g2 = scan_grid(xs2, ys2, cell_a_e, label="alpha_eps")
    g2["baseline_xy"] = (BASE["alpha"], BASE["eps"])
    grids["alpha_eps"] = g2
    plot_ab_pair(
        xs2, ys2, g2["R2_A"], g2["R2_B"], g2["kappa"],
        xlabel=r"$\alpha$", ylabel=r"$d/D$",
        title=r"$\alpha$ vs $d/D$ (fixed baseline profile)",
        baseline_xy=g2["baseline_xy"],
        out_stem="levelsets_alpha_vs_eps",
        kappa_overlay=False,  # κ fixed
    )

    # ---- 3. Amp vs φ ----
    print("\n=== Panel 3: Amp vs φ ===", flush=True)
    amp_phi = plot_amp_phi_landscape(af0, an0, mig)

    # ---- 4. D/L vs α ----
    xs4 = np.linspace(3.0, 8.0, 14)
    ys4 = np.linspace(0.5, 2.0, 12)
    print("\n=== Panel 4: D/L vs α ===", flush=True)

    def cell_dl_a(i, j, x, a):
        return eval_cell(x, BASE["f_u"], a, BASE["eps"], af_mode, mig)

    g4 = scan_grid(xs4, ys4, cell_dl_a, label="DL_alpha")
    g4["baseline_xy"] = (BASE["x"], BASE["alpha"])
    grids["DL_alpha"] = g4
    plot_ab_pair(
        xs4, ys4, g4["R2_A"], g4["R2_B"], g4["kappa"],
        xlabel=r"$D/L$", ylabel=r"$\alpha$",
        title=rf"$D/L$ vs $\alpha$ ($f_u$={BASE['f_u']}, $d/D$=0.1)",
        baseline_xy=g4["baseline_xy"],
        out_stem="levelsets_DL_vs_alpha",
        kappa_overlay=True,
    )

    # ---- summary + note ----
    print("\n=== Summary figure + NOTE ===", flush=True)
    plot_summary(grids, amp_phi)
    write_note(grids, amp_phi, verify)

    # save numeric arrays for regenerability / debugging
    summary = {
        "baseline": BASE,
        "verify": verify,
        "grids": {
            k: {
                "xs": g["xs"].tolist(),
                "ys": g["ys"].tolist(),
                "R2_A": np.asarray(g["R2_A"]).tolist(),
                "R2_B": np.asarray(g["R2_B"]).tolist(),
                "kappa": np.asarray(g["kappa"]).tolist(),
                "baseline_xy": list(g["baseline_xy"]),
                "n_cells": len(g["rows"]),
            }
            for k, g in grids.items()
        },
        "amp_phi": {
            "amps": amp_phi["amps"].tolist(),
            "phis_deg": amp_phi["phis_deg"].tolist(),
        },
        "fit_settings": {
            "n_phi_fit": N_PHI_FIT,
            "amp_grid": AMP_GRID_FIT.tolist(),
            "convention": "indep Amp_A/Amp_B shared φ (vs shared-opp); pick max min(R²)",
        },
    }
    (OUT / "levelsets_arrays.json").write_text(json.dumps(summary, indent=2))
    print("done.", flush=True)


if __name__ == "__main__":
    main()

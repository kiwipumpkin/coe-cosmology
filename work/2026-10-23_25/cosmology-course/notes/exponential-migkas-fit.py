"""Geometry-only A/A_all fit with Chris's ORIGINAL pure exponential CoE mass.

rho(r) ∝ exp(-r/L), M = M(<D); no ptshell shell, no artificial cliff.
Uses bubble_core.run (analytic Fx enclosed-mass), Ob=0.

Same Migkas Fig. 6 WLS setup as migkas-geometry-fit.py (TH1 grid, l0 shift,
chi2 with sigma=(lo+hi)/2). No sine Amp/phi — geometry only.
Free: D, LD=L/D, mc=M(<D)/M_crit(D), l0.
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
from bubble_core import IMIG, MA, MHI, MLO, RKM, GPC, TH1, run, shifted  # noqa: E402

RV = RKM / GPC  # 14.300... Gpc
OUT_DIR = Path("/workspace/cosmology-course/notes")
OUT_PNG = OUT_DIR / "exponential-migkas-fit.png"
OUT_JSON = OUT_DIR / "exponential-migkas-fit.json"
OUT_DENS = OUT_DIR / "exponential-density-bestfit.png"
OUT_CSV = OUT_DIR / "exponential-migkas-fit-curve.csv"

Y = MA
SIG = 0.5 * (MLO + MHI)
W = 1.0 / SIG ** 2
LDEG = np.loadtxt("/workspace/cosmology/migkas/fig6_digitized.csv", delimiter=",", skiprows=1, usecols=0)

# Previous ptshell WLS (geometry-only) for comparison
PTSHELL_WLS = dict(D=58.19699837268826, R2_unweighted=0.5456268309812682,
                   R2_weighted=0.5807290671992491, chi2=18.7465256529584, l0=15.1629523788697)


def r2_u(y, m):
    return float(1.0 - np.sum((y - m) ** 2) / np.sum((y - y.mean()) ** 2))


def r2_w(y, m, w):
    yb = np.sum(w * y) / np.sum(w)
    return float(1.0 - np.sum(w * (y - m) ** 2) / np.sum(w * (y - yb) ** 2))


def sky(D, LD, mc):
    """Mean-normalised a_f on 1-degree grid. None if illegal / turnaround."""
    if not (D > RV * 1.05 and 0.01 < LD < 3.0 and 1e-4 < mc < 50.0):
        return None
    DR = D / RV
    rr = run(DR, float(LD), float(mc), 0.0)  # Ob=0: pure exp, no floor
    af = rr["af"]
    if (not np.all(np.isfinite(af))) or rr["vrmin"] <= 0:
        return None
    return af / af.mean(), rr


def model_bins(D, LD, mc, l0):
    out = sky(D, LD, mc)
    if out is None:
        return None
    af, rr = out
    return shifted(af, l0), rr


def chi2(D, LD, mc, l0):
    out = model_bins(D, LD, mc, l0)
    if out is None:
        return 800.0
    pred, _ = out
    return float(np.sum(W * (pred[IMIG] - Y) ** 2))


def pack_metrics(D, LD, mc, l0, c2, rr=None):
    out = model_bins(D, LD, mc, l0)
    full, rr2 = out
    m = full[IMIG]
    if rr is None:
        rr = rr2
    L = LD * D
    return dict(
        D_Gpc=float(D),
        DR=float(D / RV),
        LD=float(LD),
        L_Gpc=float(L),
        mc=float(mc),
        Ob=0.0,
        l0=float(l0 % 360.0),
        chi2=float(c2),
        chi2_red=float(c2 / (Y.size - 4)),
        R2_unweighted=r2_u(Y, m),
        R2_weighted=r2_w(Y, m, W),
        A_min=float(full.min()),
        A_max=float(full.max()),
        vrmin=float(rr["vrmin"]),
        Ht=float(rr["Ht"]),
        contrast_Rvis=float(rr["contrast"]),  # delta at R = RV (vision)
        contrast_analytic_R3=float(2 * np.sinh(3.0 / L)),
        contrast_analytic_R14p3=float(2 * np.sinh(RV / L)),
        age_Gyr=float(rr["age"]),
    )


def objective(u):
    # u = log10 D, log10 LD, log10 mc, l0
    return chi2(10 ** u[0], 10 ** u[1], 10 ** u[2], u[3])


def fit():
    starts = [
        np.array([np.log10(30), np.log10(0.15), np.log10(0.1), 15.0]),
        np.array([np.log10(40), np.log10(0.09), np.log10(0.04), 13.0]),
        np.array([np.log10(35), np.log10(0.12), np.log10(0.08), 14.0]),
        np.array([np.log10(45), np.log10(0.08), np.log10(0.03), 15.0]),
        np.array([np.log10(25), np.log10(0.20), np.log10(0.15), 14.0]),
        np.array([np.log10(50), np.log10(0.10), np.log10(0.05), 15.2]),
        np.array([np.log10(38), np.log10(0.0877), np.log10(0.0389), 13.1]),
        np.array([np.log10(28), np.log10(0.18), np.log10(0.2), 16.0]),
    ]
    best_u, best_c = None, np.inf
    for i, u0 in enumerate(starts):
        print(f"  start {i+1}/{len(starts)}: D={10**u0[0]:.1f} LD={10**u0[1]:.4f} mc={10**u0[2]:.4g} l0={u0[3]:.1f}",
              flush=True)
        res = minimize(objective, u0, method="Nelder-Mead",
                       options=dict(maxiter=220, xatol=1e-4, fatol=1e-4, adaptive=True))
        print(f"    -> chi2={res.fun:.4f} nit={res.nit}", flush=True)
        if res.fun < best_c:
            best_c, best_u = float(res.fun), res.x.copy()
    D = 10 ** best_u[0]
    LD = 10 ** best_u[1]
    mc = 10 ** best_u[2]
    l0 = best_u[3]
    return pack_metrics(D, LD, mc, l0, best_c)


def plot_fit(fitp):
    full, _ = model_bins(fitp["D_Gpc"], fitp["LD"], fitp["mc"], fitp["l0"])
    fig, ax = plt.subplots(figsize=(10.0, 5.2))
    ax.errorbar(LDEG, Y, yerr=np.vstack([MLO, MHI]), fmt="o", ms=4.0, color="k",
                ecolor="0.45", elinewidth=0.7, capsize=0, zorder=3,
                label=r"Migkas Fig. 6  $A/A_{\mathrm{all}}$")
    ax.plot(TH1, full, color="C0", lw=2.0,
            label=(rf"WLS pure exp: $D={fitp['D_Gpc']:.1f}$ Gpc, "
                   rf"$L={fitp['L_Gpc']:.2f}$ Gpc ($L/D={fitp['LD']:.3f}$), "
                   rf"$m_c={fitp['mc']:.4g}$, $\ell_0={fitp['l0']:.1f}^\circ$"
                   f"\n"
                   rf"$R^2={fitp['R2_unweighted']:.3f}$ (unw.), "
                   rf"$\chi^2={fitp['chi2']:.2f}$; "
                   rf"ptshell WLS $R^2\approx{PTSHELL_WLS['R2_unweighted']:.2f}$, "
                   rf"$D\approx{PTSHELL_WLS['D']:.1f}$ Gpc"))
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)", fontsize=13)
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$", fontsize=13)
    ax.set_title(
        r"CoE pure exponential mass $\rho\propto e^{-r/L}$ (Ob=0; geometry only, no sine)",
        fontsize=14,
    )
    ax.legend(fontsize=14, loc="lower left", framealpha=0.92)
    ax.set_xticks(np.arange(0, 361, 45))
    ax.tick_params(labelsize=12)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=160)
    np.savetxt(OUT_CSV, np.column_stack([TH1, full]), delimiter=",",
               header="ell_deg,A_over_Aall", comments="")
    print("wrote", OUT_PNG, flush=True)


def plot_density_bestfit(fitp):
    """rho(r)/rho(D) at best-fit (D,L); M cancels. Mark Earth at D."""
    D = fitp["D_Gpc"]
    L = fitp["L_Gpc"]
    R3, R14 = 3.0, RV
    r = np.linspace(0.0, D + 1.6 * R14, 2500)
    y = np.exp(-(r - D) / L)
    fig, ax = plt.subplots(figsize=(10.0, 5.8))
    ax.plot(r, y, color="C0", lw=2.4,
            label=rf"Best-fit $\rho\propto e^{{-r/L}}$: $L={L:.2f}$ Gpc, $L/D={fitp['LD']:.3f}$")
    ax.axvline(D, color="k", lw=1.6, ls="--",
               label=rf"Earth at $D={D:.1f}$ Gpc (CoE→Earth)")
    ax.axvspan(D - R14, D + R14, color="C2", alpha=0.12,
               label=rf"Vision sphere $D\pm R$, $R={R14:.1f}$ Gpc")
    ax.axvline(D - R3, color="0.35", lw=1.0, ls=":")
    ax.axvline(D + R3, color="0.35", lw=1.0, ls=":",
               label=rf"$D\pm {R3:g}$ Gpc")
    ax.axhline(1.0, color="0.5", lw=0.8, alpha=0.7)
    d3 = fitp["contrast_analytic_R3"]
    d14 = fitp["contrast_analytic_R14p3"]
    ax.set_xlabel(r"Distance from CoE (Centre of Expansion), $r$ (Gpc)", fontsize=13)
    ax.set_ylabel(r"Mass density ratio $\rho(r)/\rho(D)$", fontsize=13)
    ax.set_title(
        rf"Best-fit pure exponential density (Migkas $A/A_{{\mathrm{{all}}}}$ WLS)"
        f"\n"
        rf"$m_c=M(<D)/M_{{\mathrm{{crit}}}}(D)={fitp['mc']:.4g}$, "
        rf"$\ell_0={fitp['l0']:.1f}^\circ$, "
        rf"$\delta(3)={d3:.3g}$, $\delta(14.3)={d14:.3g}$",
        fontsize=14,
    )
    ax.legend(fontsize=14, loc="upper right", framealpha=0.92)
    ax.set_xlim(0, D + 1.55 * R14)
    ymax = min(8.0, max(3.0, float(np.exp((D - max(D - R14, 0)) / L)) * 1.15))
    ax.set_ylim(0, ymax)
    ax.tick_params(labelsize=12)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT_DENS, dpi=160)
    print("wrote", OUT_DENS, flush=True)


def main():
    print("Fitting pure-exp geometry to Migkas Fig. 6 ...", flush=True)
    fitp = fit()
    print("BEST", {k: (round(v, 6) if isinstance(v, float) else v)
                   for k, v in fitp.items()}, flush=True)
    plot_fit(fitp)
    plot_density_bestfit(fitp)
    payload = dict(
        formula="A/A_all(ell)=a_f(ell-l0)/<a_f>, a_f=|V_source-V_Earth|/(R/D); "
                "pure rho∝exp(-r/L), Ob=0, bubble_core.run analytic Fx",
        R_Gpc=float(RV),
        profile="exp",
        free_params=["D", "LD", "mc", "l0"],
        fixed=dict(Ob=0.0, fac=1.5, taui=0.01),
        wls_fit=fitp,
        ptshell_wls_geometry_only=PTSHELL_WLS,
        note=("Amp and phi were not fitted (geometry only, same as migkas-geometry-fit). "
              "sigma=(A_err_low+A_err_high)/2. Contrast delta=2 sinh(R/L) for pure exp. "
              "mc = M(<D)/M_crit(D) today."),
    )
    OUT_JSON.write_text(json.dumps(payload, indent=1))
    print("wrote", OUT_JSON, flush=True)


if __name__ == "__main__":
    main()

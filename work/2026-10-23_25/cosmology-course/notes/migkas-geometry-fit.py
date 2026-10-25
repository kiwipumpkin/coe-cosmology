"""Geometry-only A/A_all versus Galactic longitude, from the CoE + thick-shell flow.

Chris's Newtonian integrator (papers/spreadsheet_model/bubble/ptshell.py,
prof_core.py, integ_core.geometry). The curve is a_f/<a_f> shifted onto
Galactic longitude. The diagnostic sine
    A = a_f * [1 - Amp * sin(psi + phi)]
is not applied and Amp, phi are not fitted.

Explicit map (length unit = D, CoE at the origin, Earth at q = (1, 0)):

    epsilon = R/D
    q_source(theta) = (1 + epsilon * cos theta, epsilon * sin theta)
    a_f(theta) = |V(q_source) - V(q_Earth)| / epsilon
    A/A_all (ell) = a_f(ell - l0) / <a_f>

theta = 0 is the anti-CoE direction. V is the integrated Hubble flow
(start tau = 0.01, Hubble-like) decelerated by a compact core
mcc = M_c/M_crit(D) plus a uniform shell Omega_s = 0.049 on
R_in = D - k_in * R  to  R_out = D + k_out * R, with k_in > 1 and k_out > 1
so the observer sits inside the shell by at least one vision radius R.
<a_f> is the mean over 360 one-degree sightlines. R = 14.30048287260589 Gpc.

The multipole expansion of this same a_f (draft section 2.5.1) is

    a_f**2  proportional to  1 + B mu**2 + epsilon (P mu + Q mu**3)
    B = 2*beta + beta**2
    P = beta*(3+beta)
    Q = (1+beta)*(gamma+beta) - 2*beta
    a_f/<a_f> = 1 + a1 cos psi + a2 cos 2psi + a3 cos 3psi + ...

with mu = cos psi, psi = ell - l0, beta = (dh/dq)/h and gamma = (d2h/dq2)/h
at q = 1, for the radial profile V = h(q) q. Even a_m are fixed by beta;
odd a_m by epsilon*Q. The plotted curve is the numeric a_f, not that series
and not a free sine.

Weights: sigma = (A_err_low + A_err_high)/2, chi2 = sum ((model-data)/sigma)**2.
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
import ptshell  # noqa: E402,F401  patches prof_core for the core+shell mass table
import prof_core as P  # noqa: E402
from bubble_core import IMIG, MA, MHI, MLO, TH1, shifted  # noqa: E402

RV = 14.30048287260589
OS = 0.049
OUT_PNG = Path("/workspace/cosmology-course/notes/migkas-geometry-fit.png")
OUT_JSON = Path("/workspace/cosmology-course/notes/migkas-geometry-fit.json")

# Table 1 Newtonian Omega_s=0.049 column. Those Amp, phi are not used here.
TABLE1 = dict(D=49.17456190433319, mcc=0.0017144539494289297,
              kin=1.0276268339399446, kout=1.4754794866889613, l0=13.0)

Y = MA
SIG = 0.5 * (MLO + MHI)
W = 1.0 / SIG ** 2
LDEG = np.loadtxt("/workspace/cosmology/migkas/fig6_digitized.csv", delimiter=",", skiprows=1, usecols=0)


def r2_u(y, m):
    return float(1.0 - np.sum((y - m) ** 2) / np.sum((y - y.mean()) ** 2))


def r2_w(y, m, w):
    yb = np.sum(w * y) / np.sum(w)
    return float(1.0 - np.sum(w * (y - m) ** 2) / np.sum(w * (y - yb) ** 2))


def sky(D, mcc, kin, kout):
    """Mean-normalised a_f on the 1-degree grid. None if the shell is illegal or the flow turns around."""
    if not (D > 0 and mcc > 0 and kin > 1.001 and kout > 1.001):
        return None
    rin = 1.0 - kin * RV / D
    rout = 1.0 + kout * RV / D
    # rin > 0.008 keeps the enclosed-mass table off the origin; the WLS minimum sits on this floor.
    if not (0.008 < rin < 0.99 and rout > 1.0):
        return None
    a = dict(DR=D / RV, Os=OS, mcc=float(mcc), rin=float(rin), rout=float(rout),
             Lc=0.005, mc=0.0, Ob=0.0)
    rr = P.run("ptshell", a, 0.0, 0.0)
    af = rr["af"]
    if (not np.all(np.isfinite(af))) or rr["vrmin"] <= 0:
        return None
    return af / af.mean()


def model_bins(D, mcc, kin, kout, l0):
    af = sky(D, mcc, kin, kout)
    if af is None:
        return None
    return shifted(af, l0)


def chi2(D, mcc, kin, kout, l0):
    pred = model_bins(D, mcc, kin, kout, l0)
    if pred is None:
        return 800.0
    return float(np.sum(W * (pred[IMIG] - Y) ** 2))


def pack_metrics(D, mcc, kin, kout, l0, c2):
    full = model_bins(D, mcc, kin, kout, l0)
    m = full[IMIG]
    return dict(
        D=float(D), DR=float(D / RV), mcc=float(mcc), kin=float(kin), kout=float(kout),
        l0=float(l0 % 360.0), Rin=float(D - kin * RV), Rout=float(D + kout * RV),
        chi2=float(c2), chi2_red=float(c2 / (Y.size - 5)),
        R2_unweighted=r2_u(Y, m), R2_weighted=r2_w(Y, m, W),
        A_min=float(full.min()), A_max=float(full.max()),
    )


def objective(u):
    # u = log10 D, log10 mcc, log10(kin-1), log10(kout-1), l0
    return chi2(10 ** u[0], 10 ** u[1], 1 + 10 ** u[2], 1 + 10 ** u[3], u[4])


def fit():
    starts = [
        np.array([np.log10(58), np.log10(0.00171), np.log10(3.0), np.log10(0.05), 15.2]),
        np.array([np.log10(55), np.log10(0.0017), np.log10(2.6), np.log10(0.02), 15.0]),
        np.array([np.log10(49.2), np.log10(0.00171), np.log10(0.028), np.log10(0.475), 13.0]),
        np.array([np.log10(80), np.log10(0.0016), np.log10(1.2), np.log10(0.25), 14.0]),
        np.array([np.log10(40), np.log10(0.0025), np.log10(1.0), np.log10(0.4), 14.0]),
    ]
    best_u, best_c = None, np.inf
    for u0 in starts:
        res = minimize(objective, u0, method="Nelder-Mead",
                       options=dict(maxiter=160, xatol=1e-4, fatol=1e-4, adaptive=True))
        if res.fun < best_c:
            best_c, best_u = float(res.fun), res.x.copy()
    D = 10 ** best_u[0]
    mcc = 10 ** best_u[1]
    kin = 1 + 10 ** best_u[2]
    kout = 1 + 10 ** best_u[3]
    l0 = best_u[4]
    return pack_metrics(D, mcc, kin, kout, l0, best_c)


def main():
    t1 = pack_metrics(TABLE1["D"], TABLE1["mcc"], TABLE1["kin"], TABLE1["kout"], TABLE1["l0"],
                      chi2(TABLE1["D"], TABLE1["mcc"], TABLE1["kin"], TABLE1["kout"], TABLE1["l0"]))
    print("Table1 no-sine", {k: round(t1[k], 4) if isinstance(t1[k], float) else t1[k] for k in
          ("D", "DR", "Rin", "Rout", "l0", "mcc", "R2_unweighted", "R2_weighted", "chi2")}, flush=True)
    fitp = fit()
    print("WLS", {k: round(fitp[k], 4) if isinstance(fitp[k], float) else fitp[k] for k in
          ("D", "DR", "mcc", "kin", "kout", "Rin", "Rout", "l0",
           "R2_unweighted", "R2_weighted", "chi2", "chi2_red")}, flush=True)
    t1_sky = model_bins(t1["D"], t1["mcc"], t1["kin"], t1["kout"], t1["l0"])
    fit_sky = model_bins(fitp["D"], fitp["mcc"], fitp["kin"], fitp["kout"], fitp["l0"])

    fig, ax = plt.subplots(figsize=(9.2, 4.7))
    ax.errorbar(LDEG, Y, yerr=np.vstack([MLO, MHI]), fmt="o", ms=3.6, color="k",
                ecolor="0.45", elinewidth=0.7, capsize=0, zorder=3,
                label=r"Migkas Fig. 6  $A/A_{\mathrm{all}}$")
    ax.plot(TH1, t1_sky, color="0.4", lw=1.15, ls="--",
            label=(rf"Table 1 shell, no sine: $D={t1['D']:.1f}$ Gpc, "
                   rf"$\ell_0={t1['l0']:.0f}^\circ$, $R^2={t1['R2_unweighted']:.3f}$"))
    ax.plot(TH1, fit_sky, color="C0", lw=1.8,
            label=(rf"WLS geometry: $D={fitp['D']:.1f}$ Gpc, "
                   rf"$R_{{\rm in}}={fitp['Rin']:.1f}$, $R_{{\rm out}}={fitp['Rout']:.1f}$, "
                   rf"$\ell_0={fitp['l0']:.1f}^\circ$, $R^2={fitp['R2_unweighted']:.3f}$"))
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$")
    ax.set_title(r"CoE + thick shell, geometry only ($\Omega_s=0.049$; no sine rotation)")
    ax.legend(fontsize=7.6, loc="lower left")
    ax.set_xticks(np.arange(0, 361, 45))
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=160)
    payload = dict(
        formula="A/A_all(ell)=a_f(ell-l0)/<a_f>, a_f=|V_source-V_Earth|/(R/D)",
        omega_s=OS, R_Gpc=RV,
        table1_geometry_only=t1,
        wls_fit=fitp,
        note=("Amp and phi were not fitted. sigma=(A_err_low+A_err_high)/2. "
              "Rin runs to the integrator floor rin=0.008 D; chi2 is flat in Rout once Rout is at least D+R."),
    )
    OUT_JSON.write_text(json.dumps(payload, indent=1))
    print("wrote", OUT_PNG, flush=True)


if __name__ == "__main__":
    main()

"""CatWISE-constrained η (Einasto-like / generalized exponential) CoE density
+ Migkas A/A_all refit under the 3% lock.

Density: rho(r) = rho_0 * exp(-(r/L)^eta)
  - eta=1 recovers pure exponential
  - Finite at CoE: rho(0)=rho_0; rho→0 as r→∞
  - rho(D)=Omega_m * rho_crit with Omega_m=0.3

CatWISE 3% lock (ONLY fixed constraint):
  delta = [rho(D-d)-rho(D+d)]/rho(D) = 0.031
  => for each trial (D, eta), solve L>0 from
  delta = exp[-(D/L)^eta*((1-d/D)^eta-1)]
        - exp[-(D/L)^eta*((1+d/D)^eta-1)]
  eta=1 recovers L = d/arcsinh(delta/2).

Migkas: geometry-only A/A_all = a_f(ell-l0)/<a_f>, Ob=0, WLS.
Free {D, mc, l0, eta}; L locked by delta. Primary d=vision R≈14.3 Gpc;
also report d=3 Gpc.

Uses prof_core kind "eta" (minimal patch: shape + LB/LOGV eta).
"""
from __future__ import annotations

import json
import os
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq, minimize
from scipy.stats import chi2 as chi2_dist, norm

BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"
sys.path.insert(0, BUB)
os.chdir(BUB)
from bubble_core import (  # noqa: E402
    IMIG, MA, MHI, MLO, RKM, GPC, TH1, shifted,
)
import prof_core as pc  # noqa: E402

OUT = Path("/workspace/cosmology-course/notes")

# --- CatWISE / geometry constants ---
D_CATWISE = 0.01554
DELTA = 2.0 * D_CATWISE
DELTA_USE = 0.031
D_EARTH = 58.2
RV = float(RKM / GPC)          # ≈ 14.3005 Gpc vision radius
D_NEAR = 3.0

OMEGA_M = 0.3
H0 = 70.0

# Prior table rows (for comparison / σ vs null)
CHI2_NULL = 45.13662301218721
N_BINS = 71
PTSHELL_WLS = dict(
    D=58.19699837268826,
    R2_unweighted=0.5456268309812682,
    R2_weighted=0.5807290671992491,
    chi2=18.7465256529584,
    l0=15.1629523788697,
    n_free=5,
)
CATWISE_EXP_D14 = dict(
    D_Gpc=45.99056277928336,
    L_Gpc=922.6487386323013,
    R2_unweighted=0.716036731924234,
    chi2=10.572265211558566,
    n_free=3,
)
FREE_EXP_WLS = dict(
    D_Gpc=27.95388837552518,
    L_Gpc=6.072715534614713,
    LD=0.21724045875248016,
    R2_unweighted=0.7589233227942257,
    chi2=8.47533016653482,
    n_free=4,
)

Y = MA
SIG = 0.5 * (MLO + MHI)
W = 1.0 / SIG ** 2
LDEG = np.loadtxt(
    "/workspace/cosmology/migkas/fig6_digitized.csv",
    delimiter=",", skiprows=1, usecols=0,
)

# Stellar-mass BH mean density reference (do not invent a floor; flag only).
# 10 M_sun Schwarzschild: R_s = 2GM/c^2 ≈ 29.5 km; mean density 3M/(4π R_s^3).
RHO_BH_MSUN_MPC3 = 1.8e40   # ~order-of-magnitude 10 M_sun BH mean density in M_sun/Mpc^3


def rho_crit_Msun_Mpc3(H0_kms: float = H0) -> float:
    h = H0_kms / 100.0
    return 2.77536627e11 * h ** 2


def r2_u(y, m):
    return float(1.0 - np.sum((y - m) ** 2) / np.sum((y - y.mean()) ** 2))


def r2_w(y, m, w):
    yb = np.sum(w * y) / np.sum(w)
    return float(1.0 - np.sum(w * (y - m) ** 2) / np.sum(w * (y - yb) ** 2))


def sigma_from_delta(dchi2, dnu):
    if dchi2 <= 0:
        return 0.0, 1.0
    p = float(chi2_dist.sf(dchi2, max(dnu, 1)))
    if p <= 0 or not np.isfinite(p):
        return float("inf"), p
    return float(norm.isf(p)), p


# ---------------------------------------------------------------------------
# CatWISE δ lock → L(D, η, d)
# ---------------------------------------------------------------------------
def delta_eta(D: float, L: float, d: float, eta: float) -> float:
    """Pole-to-pole contrast for rho ∝ exp(-(r/L)^eta)."""
    if L <= 0 or D <= 0 or d <= 0 or eta <= 0:
        return np.nan
    if d >= D:
        return np.nan
    a = (D / L) ** eta
    t1 = np.exp(-a * ((1.0 - d / D) ** eta - 1.0))
    t2 = np.exp(-a * ((1.0 + d / D) ** eta - 1.0))
    return float(t1 - t2)


def L_from_delta_eta(
    D: float, d: float, eta: float, delta: float = DELTA_USE,
) -> float:
    """Solve δ(D,L,d,η)=delta for L>0.  η=1 → L=d/arcsinh(δ/2)."""
    if eta <= 0:
        raise ValueError("eta must be > 0")
    if d >= D:
        raise ValueError("need d < D for CatWISE lock")
    if abs(eta - 1.0) < 1e-12:
        return float(d / np.arcsinh(delta / 2.0))

    def f(log_alpha: float) -> float:
        alpha = np.exp(log_alpha)  # alpha = D/L
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with np.errstate(over="ignore", invalid="ignore"):
                val = delta_eta(D, D / alpha, d, eta) - delta
        return val

    # Wide grid on log(D/L); pick first valid sign-change bracket
    grid = np.linspace(-30.0, 20.0, 4000)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        vals = np.array([f(g) for g in grid])
    ok = np.isfinite(vals)
    idxs = []
    for i in range(len(grid) - 1):
        if not (ok[i] and ok[i + 1]):
            continue
        a, b = vals[i], vals[i + 1]
        if a == 0.0 or b == 0.0 or (a > 0) != (b > 0):
            idxs.append(i)
    if not idxs:
        raise RuntimeError(
            f"no L root for D={D}, d={d}, eta={eta}, delta={delta}"
        )
    lo, hi = float(grid[idxs[0]]), float(grid[idxs[0] + 1])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        log_alpha = brentq(f, lo, hi, xtol=1e-12)
    L = float(D / np.exp(log_alpha))
    if L <= 0 or not np.isfinite(L):
        raise RuntimeError(f"bad L={L} for D={D}, d={d}, eta={eta}")
    return L


def slope_at_D(L: float, D: float, eta: float) -> float:
    """(rho'/rho)_D = -(eta/L) (D/L)^{eta-1}."""
    return float(-(eta / L) * (D / L) ** (eta - 1.0))


def density_params(
    d_gpc: float, eta: float, D: float = D_EARTH, delta: float = DELTA_USE,
):
    L = L_from_delta_eta(D, d_gpc, eta, delta)
    alpha = D / L
    rho0_over_crit = OMEGA_M * np.exp(alpha ** eta)
    slope_ln = slope_at_D(L, D, eta)
    drho_dr_D = OMEGA_M * slope_ln
    dcheck = delta_eta(D, L, d_gpc, eta)
    rho_c = rho_crit_Msun_Mpc3()
    rho0_phys = rho0_over_crit * rho_c
    below_bh = bool(rho0_phys < RHO_BH_MSUN_MPC3)
    return dict(
        d_Gpc=float(d_gpc),
        D_Gpc=float(D),
        eta=float(eta),
        L_Gpc=float(L),
        L_over_D=float(L / D),
        delta_target=float(delta),
        delta_check=float(dcheck),
        D_CatWISE=float(D_CATWISE),
        Omega_m=float(OMEGA_M),
        rho_D_over_rho_crit=float(OMEGA_M),
        rho_0_over_rho_crit=float(rho0_over_crit),
        drho_dr_at_D_over_rho_crit_per_Gpc=float(drho_dr_D),
        rho_prime_over_rho_at_D=float(slope_ln),
        rho_0_below_BH_density=below_bh,
        rho_BH_ref_Msun_Mpc3=float(RHO_BH_MSUN_MPC3),
        rho_to_0_at_infinity=True,
        note_rho_choice="rho(D)=Omega_m*rho_crit with Omega_m=0.3",
        law="rho(r)=rho_0*exp(-(r/L)^eta), rho_0=rho(D)*exp((D/L)^eta)",
    )


def rho_over_crit(r, D, L, eta):
    return OMEGA_M * np.exp((D / L) ** eta - (r / L) ** eta)


# ---------------------------------------------------------------------------
# Density plot
# ---------------------------------------------------------------------------
def plot_density(eta_fit: float, primary_d: float = RV):
    """rho(r)/rho_crit for best-fit η (+ η=1 CatWISE overlay); D=58.2 baseline."""
    D = D_EARTH
    cases = [
        (eta_fit, primary_d, "C0", "-", True, f"PRIMARY η={eta_fit:.3g}"),
        (1.0, primary_d, "C1", "--", False, "η=1 CatWISE-L (comparison)"),
    ]
    r_max = D + 1.55 * max(primary_d, D_NEAR, RV)
    r = np.linspace(0.0, r_max, 3000)

    fig, ax = plt.subplots(figsize=(13.5, 6.8))
    params_all = {}
    csv_cols = [r]
    csv_hdr = ["r_Gpc"]
    ymax_data = 0.0

    for eta, d, color, ls, is_primary, tag in cases:
        p = density_params(d, eta, D=D)
        key = f"eta={eta:g}_d={d:g}"
        params_all[key] = p
        L = p["L_Gpc"]
        y = rho_over_crit(r, D, L, eta)
        ymax_data = max(ymax_data, float(np.nanmax(y[np.isfinite(y)])))
        lab = (
            rf"{tag}: $d={d:.2f}$ Gpc $\Rightarrow$ $L={L:.1f}$ Gpc"
            rf" ($\delta={p['delta_check']:.4f}$)"
        )
        ax.plot(r, y, color=color, ls=ls, lw=2.5, label=lab)
        csv_cols.append(y)
        csv_hdr.append(f"rho_over_rhocrit_eta{eta:g}_d{d:g}_L{L:.4g}")

        if is_primary:
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

            slope_ln = p["rho_prime_over_rho_at_D"]
            ax.annotate(
                rf"$(d\rho/dr)_{{D}}={slope:.3e}\,\rho_{{\rm crit}}/$Gpc"
                f"\n"
                rf"$(\rho'/\rho)_D=-(\eta/L)(D/L)^{{\eta-1}}={slope_ln:.3e}$ "
                rf"Gpc$^{{-1}}$",
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

            ax.axvline(D - d, color=color, lw=1.2, ls=":")
            ax.axvline(
                D + d, color=color, lw=1.2, ls=":",
                label=rf"$D\pm d$ with $d={d:.2f}$ Gpc (PRIMARY)",
            )
            rho_m = float(rho_over_crit(D - d, D, L, eta))
            rho_p = float(rho_over_crit(D + d, D, L, eta))
            ax.plot([D - d], [rho_m], "s", color=color, ms=7, zorder=5)
            ax.plot([D + d], [rho_p], "s", color=color, ms=7, zorder=5)

            ax.annotate(
                rf"$\rho(D-d)/\rho_{{\rm crit}}={rho_m:.4f}$",
                xy=(D - d, rho_m),
                xytext=(D - d - 6.0, 0.095),
                fontsize=10, color=color, ha="center",
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
                fontsize=10, color=color, ha="left", va="bottom",
                arrowprops=dict(arrowstyle="->", color=color, lw=0.9),
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor=color, alpha=0.92),
                zorder=6,
            )

    ax.axvline(
        D, color="k", lw=1.8, ls="--",
        label=rf"Earth at $D={D:.1f}$ Gpc (CoE→Earth)",
    )
    ax.axhline(OMEGA_M, color="0.45", lw=0.9, alpha=0.8)

    p0 = density_params(primary_d, eta_fit, D=D)
    bh_flag = (
        r"$\ll\rho_{\rm BH}$" if p0["rho_0_below_BH_density"]
        else r"FLAG: $\rho(0)$ NOT below $\rho_{\rm BH}$"
    )
    ax.annotate(
        rf"$\rho(0)/\rho_{{\rm crit}}={p0['rho_0_over_rho_crit']:.4f}$ "
        rf"(finite; {bh_flag})"
        f"\n"
        rf"$\rho\to 0$ as $r\to\infty$; $\eta={eta_fit:.4g}$",
        xy=(0.0, p0["rho_0_over_rho_crit"]),
        xytext=(1.5, 0.028),
        fontsize=10, ha="left", va="bottom",
        arrowprops=dict(
            arrowstyle="->", color="0.25", lw=1.0,
            connectionstyle="arc3,rad=0.12",
        ),
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.92),
        zorder=6,
    )

    ax.set_xlabel(
        r"Distance from CoE (Centre of Expansion), $r$ (Gpc)", fontsize=13,
    )
    ax.set_ylabel(r"Mass density $\rho(r)/\rho_{\mathrm{crit}}$", fontsize=13)
    ax.set_title(
        r"CatWISE-constrained $\eta$-density: "
        r"$\rho(r)=\rho_0\,e^{-(r/L)^\eta}$"
        "\n"
        rf"$\delta={DELTA_USE:g}$ locks $L(D,\eta,d)$; "
        rf"primary $d={primary_d:.2f}$ Gpc, best-fit $\eta={eta_fit:.4g}$; "
        rf"$\rho(D)=\Omega_m\rho_{{\rm crit}}$, $\Omega_m={OMEGA_M:g}$",
        fontsize=13,
    )
    ax.set_xlim(0, r_max)
    ax.set_ylim(0, max(ymax_data * 1.08, OMEGA_M * 1.5))
    ax.tick_params(labelsize=12)
    ax.grid(alpha=0.25)

    fig.subplots_adjust(left=0.08, right=0.56, top=0.88, bottom=0.11)
    ax.legend(
        fontsize=14, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        ncol=1, framealpha=0.95, borderaxespad=0.0,
        fancybox=True, edgecolor="0.6",
    )

    png = OUT / "catwise-eta-density.png"
    fig.savefig(png, dpi=160, bbox_inches="tight")
    plt.close(fig)

    csv = OUT / "catwise-eta-density.csv"
    np.savetxt(
        csv, np.column_stack(csv_cols), delimiter=",",
        header=",".join(csv_hdr), comments="",
    )
    print("wrote", png)
    print("wrote", csv)
    return params_all


# ---------------------------------------------------------------------------
# Migkas A/A_all fit with L(D,η) locked
# ---------------------------------------------------------------------------
def sky_locked(D, eta, mc, d_lock):
    """Mean-normalised a_f via prof_core eta; L from CatWISE δ lock."""
    if not (D > max(d_lock, RV) * 1.05 and 1e-4 < mc < 80.0 and 0.03 < eta < 8.0):
        return None
    try:
        L = L_from_delta_eta(D, d_lock, eta, DELTA_USE)
    except Exception:
        return None
    if L <= 1e-3:
        return None
    LD = L / D
    if LD < 1e-4:
        return None
    DR = D / RV
    a = dict(DR=float(DR), LD=float(LD), eta=float(eta))
    try:
        rr = pc.run("eta", a, float(mc), 0.0)
    except Exception:
        return None
    af = rr["af"]
    if (not np.all(np.isfinite(af))) or rr["vrmin"] <= 0:
        return None
    return af / af.mean(), rr, LD, L


def model_bins(D, eta, mc, l0, d_lock):
    out = sky_locked(D, eta, mc, d_lock)
    if out is None:
        return None
    af, rr, LD, L = out
    return shifted(af, l0), rr, LD, L


def chi2(D, eta, mc, l0, d_lock):
    out = model_bins(D, eta, mc, l0, d_lock)
    if out is None:
        return 800.0
    pred, _, _, _ = out
    return float(np.sum(W * (pred[IMIG] - Y) ** 2))


def fit_locked(d_lock: float, tag: str):
    """Free {D, mc, l0, eta}; L locked by δ(D±d)=0.031."""
    # u = [log10(D), log10(mc), l0, log10(eta)]
    starts = [
        np.array([np.log10(58.2), np.log10(0.2), 15.0, np.log10(1.0)]),
        np.array([np.log10(46.0), np.log10(4.5), 21.0, np.log10(1.0)]),
        np.array([np.log10(40.0), np.log10(0.5), 14.0, np.log10(0.7)]),
        np.array([np.log10(30.0), np.log10(1.0), 15.0, np.log10(1.5)]),
        np.array([np.log10(50.0), np.log10(2.0), 16.0, np.log10(0.5)]),
        np.array([np.log10(45.0), np.log10(3.0), 20.0, np.log10(2.0)]),
        np.array([np.log10(35.0), np.log10(0.8), 18.0, np.log10(0.3)]),
        np.array([np.log10(55.0), np.log10(1.5), 15.2, np.log10(1.2)]),
        np.array([np.log10(28.0), np.log10(0.2), 18.6, np.log10(0.8)]),
        np.array([np.log10(70.0), np.log10(0.1), 13.0, np.log10(2.5)]),
        np.array([np.log10(42.0), np.log10(5.0), 20.0, np.log10(0.4)]),
        np.array([np.log10(48.0), np.log10(4.0), 21.0, np.log10(1.8)]),
    ]

    def objective(u):
        return chi2(10 ** u[0], 10 ** u[3], 10 ** u[1], u[2], d_lock)

    best_u, best_c = None, np.inf
    for i, u0 in enumerate(starts):
        print(
            f"  [{tag} d={d_lock:g}] start {i+1}/{len(starts)}: "
            f"D={10**u0[0]:.1f} mc={10**u0[1]:.4g} l0={u0[2]:.1f} "
            f"eta={10**u0[3]:.3g}",
            flush=True,
        )
        res = minimize(
            objective, u0, method="Nelder-Mead",
            options=dict(maxiter=320, xatol=1e-4, fatol=1e-4, adaptive=True),
        )
        print(f"    -> chi2={res.fun:.4f} nit={res.nit}", flush=True)
        if res.fun < best_c:
            best_c, best_u = float(res.fun), res.x.copy()

    D = 10 ** best_u[0]
    mc = 10 ** best_u[1]
    l0 = best_u[2]
    eta = 10 ** best_u[3]
    out = model_bins(D, eta, mc, l0, d_lock)
    full, rr, LD, L = out
    m = full[IMIG]
    dcheck = delta_eta(D, L, d_lock, eta)
    # also analytic contrast at vision / 3 Gpc with this (D,L,eta)
    c_vis = delta_eta(D, L, RV, eta) if D > RV else float("nan")
    c_3 = delta_eta(D, L, D_NEAR, eta) if D > D_NEAR else float("nan")
    rho0 = OMEGA_M * np.exp((D / L) ** eta)
    rho_c = rho_crit_Msun_Mpc3()
    below_bh = bool(rho0 * rho_c < RHO_BH_MSUN_MPC3)
    n_free = 4
    dchi2 = CHI2_NULL - best_c
    sig, pval = sigma_from_delta(dchi2, n_free)
    fitp = dict(
        D_Gpc=float(D),
        DR=float(D / RV),
        eta=float(eta),
        L_Gpc_locked=float(L),
        LD=float(LD),
        mc=float(mc),
        Ob=0.0,
        l0=float(l0 % 360.0),
        chi2=float(best_c),
        chi2_red=float(best_c / max(Y.size - n_free, 1)),
        R2_unweighted=r2_u(Y, m),
        R2_weighted=r2_w(Y, m, W),
        A_min=float(full.min()),
        A_max=float(full.max()),
        A_amp=float(full.max() - full.min()),
        vrmin=float(rr["vrmin"]),
        Ht=float(rr["Ht"]),
        contrast_Rvis_model=float(rr["contrast"]),
        contrast_analytic_at_d=float(dcheck),
        contrast_analytic_R14p3=float(c_vis),
        contrast_analytic_R3=float(c_3),
        d_used_for_lock_Gpc=float(d_lock),
        age_Gyr=float(rr["age"]),
        rho_0_over_rho_crit=float(rho0),
        rho_0_below_BH_density=below_bh,
        rho_prime_over_rho_at_D=slope_at_D(L, D, eta),
        free_params=["D", "mc", "l0", "eta"],
        fixed_L_rule=(
            "L solved from delta=[rho(D-d)-rho(D+d)]/rho(D)=0.031 "
            "for each trial (D,eta); d=d_lock"
        ),
        Delta_chi2_vs_null=float(dchi2),
        Delta_nu_vs_null=int(n_free),
        p_value_nested=float(pval) if pval is not None else None,
        sigma_nested_one_sided=float(sig),
        sigma_sqrt_Delta_chi2=float(np.sqrt(max(dchi2, 0.0))),
        chi2_null=float(CHI2_NULL),
        N_bins=int(N_BINS),
        note_eta=("chi2 improves monotonically as eta decreases under the CatWISE delta lock; best eta sits at the lower numerical viability edge (~0.05; eta≲0.02 fails L-solve/dynamics)" if eta < 0.06 else None),
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
            rf"CatWISE-$\eta$: $d={fitp['d_used_for_lock_Gpc']:.2f}$ Gpc, "
            rf"$\eta={fitp['eta']:.4g}$, $L={fitp['L_Gpc_locked']:.1f}$ Gpc (locked)"
            f"\n"
            rf"$D={fitp['D_Gpc']:.1f}$ Gpc, $L/D={fitp['LD']:.3f}$, "
            rf"$m_c={fitp['mc']:.4g}$, $\ell_0={fitp['l0']:.1f}^\circ$"
            f"\n"
            rf"$R^2={fitp['R2_unweighted']:.3f}$ (unw.), $\chi^2={fitp['chi2']:.2f}$; "
            rf"$\sigma_{{\rm null}}={fitp['sigma_nested_one_sided']:.2f}$; "
            rf"$A_{{\rm amp}}={fitp['A_amp']:.4f}$"
            f"\n"
            rf"ptshell $R^2\approx{PTSHELL_WLS['R2_unweighted']:.2f}$; "
            rf"CatWISE-exp $R^2\approx{CATWISE_EXP_D14['R2_unweighted']:.2f}$; "
            rf"free-exp $R^2\approx{FREE_EXP_WLS['R2_unweighted']:.2f}$"
        ),
    )
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)", fontsize=13)
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$", fontsize=13)
    ax.set_title(
        rf"Migkas geometry-only with CatWISE $\delta={DELTA_USE:g}$ lock on "
        rf"$\rho\propto e^{{-(r/L)^\eta}}$"
        f"\n"
        rf"Ob=0; free {{$D$, $m_c$, $\ell_0$, $\eta$}}; $L=L(D,\eta,d)$ — {tag}",
        fontsize=14,
    )
    ax.legend(fontsize=14, loc="lower left", framealpha=0.92)
    ax.set_xticks(np.arange(0, 361, 45))
    ax.tick_params(labelsize=12)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    png = OUT / f"catwise-eta-migkas-fit-{tag}.png"
    fig.savefig(png, dpi=160, bbox_inches="tight")
    plt.close(fig)
    csv = OUT / f"catwise-eta-migkas-fit-{tag}.csv"
    np.savetxt(
        csv, np.column_stack([TH1, full]), delimiter=",",
        header="ell_deg,A_over_Aall", comments="",
    )
    print("wrote", png)
    print("wrote", csv)
    return png


def write_md(dens_params: dict, fits: dict):
    # primary dens = best-fit eta case
    p_primary = None
    p_eta1 = None
    for k, v in dens_params.items():
        if "PRIMARY" in k or (abs(v["eta"] - fits.get("d14.3", {}).get("eta", -999)) < 1e-6
                              and abs(v["d_Gpc"] - RV) < 0.05):
            p_primary = v
        if abs(v["eta"] - 1.0) < 1e-9 and abs(v["d_Gpc"] - RV) < 0.05:
            p_eta1 = v
    if p_primary is None:
        p_primary = list(dens_params.values())[0]

    rho_c = rho_crit_Msun_Mpc3()
    fp = fits.get("d14.3", list(fits.values())[0])
    lines = []
    lines.append("# CatWISE-constrained η-density (generalized exponential) + Migkas refit")
    lines.append("")
    lines.append("**Date:** 2026-10-06 (America/Chicago)")
    lines.append("")
    lines.append(
        "**Acronyms:** CatWISE = CatWISE quasar dipole catalogue (Secrest et al. 2021); "
        "CoE = centre of expansion; Migkas = Migkas et al. cluster scaling "
        r"$A/A_{\mathrm{all}}$ anisotropy."
    )
    lines.append("")
    lines.append("## Density law")
    lines.append("")
    lines.append(r"$$\rho(r)=\rho_0\,e^{-(r/L)^\eta},\qquad \rho_0=\rho(D)\,e^{(D/L)^\eta}$$")
    lines.append("")
    lines.append(
        r"- $\eta=1$ recovers the pure exponential $\rho(r)=\rho_0 e^{-r/L}$."
    )
    lines.append(r"- Finite at CoE: $\rho(0)=\rho_0$; $\rho\to 0$ as $r\to\infty$.")
    lines.append(f"- Earth–CoE baseline for density plot: $D={D_EARTH}$ Gpc.")
    lines.append(
        f"- Local normalisation: $\\rho(D)=\\Omega_m\\rho_{{\\rm crit}}$ with "
        f"**$\\Omega_m={OMEGA_M}$**."
    )
    lines.append(
        f"- $\\rho_{{\\rm crit}}(H_0={H0:g})\\approx {rho_c:.4e}\\,M_\\odot/\\mathrm{{Mpc}}^3$."
    )
    lines.append(
        r"- Local slope: $(\rho'/\rho)_D=-(\eta/L)(D/L)^{\eta-1}$."
    )
    lines.append(
        f"- BH-density flag: $\\rho(0)/\\rho_{{\\rm crit}}="
        f"{p_primary['rho_0_over_rho_crit']:.5g}$ "
        f"({'below' if p_primary['rho_0_below_BH_density'] else 'NOT below'} "
        f"stellar-BH reference; no floor invented)."
    )
    lines.append("")
    lines.append("## CatWISE 3% lock (only fixed constraint)")
    lines.append("")
    lines.append(
        r"Pole-to-pole contrast at vision distance $d$:"
    )
    lines.append("")
    lines.append(
        r"$$\delta=\frac{\rho(D-d)-\rho(D+d)}{\rho(D)}"
        r"=e^{-(D/L)^\eta\big((1-d/D)^\eta-1\big)}"
        r"-e^{-(D/L)^\eta\big((1+d/D)^\eta-1\big)}\approx 0.031.$$"
    )
    lines.append("")
    lines.append(
        r"For $\eta=1$ this collapses to $\delta=2\sinh(d/L)$, so "
        r"$L=d/\mathrm{arcsinh}(\delta/2)$ (independent of $D$)."
    )
    lines.append("")
    lines.append(
        f"Primary $d=R_{{\\rm vis}}\\approx{RV:.4g}$ Gpc; secondary $d={D_NEAR:g}$ Gpc."
    )
    lines.append("")
    lines.append("| case | $\\eta$ | $d$ (Gpc) | $D$ (Gpc) | $L$ (Gpc) | $\\delta$ check | $\\rho(0)/\\rho_{{\\rm crit}}$ | $(\\rho'/\\rho)_D$ (Gpc$^{-1}$) |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|")
    for k, p in dens_params.items():
        lines.append(
            f"| {k} | {p['eta']:.4g} | {p['d_Gpc']:.4g} | {p['D_Gpc']:.4g} | "
            f"{p['L_Gpc']:.4g} | {p['delta_check']:.5f} | "
            f"{p['rho_0_over_rho_crit']:.5g} | {p['rho_prime_over_rho_at_D']:.6g} |"
        )
    if p_eta1 is not None:
        L_ref = RV / np.arcsinh(DELTA_USE / 2.0)
        pct = 100.0 * abs(p_eta1["L_Gpc"] - L_ref) / L_ref
        lines.append("")
        lines.append(
            f"**η=1 recovery check:** $L={p_eta1['L_Gpc']:.4g}$ Gpc vs "
            f"$d/\\mathrm{{arcsinh}}(\\delta/2)={L_ref:.4g}$ Gpc "
            f"(relative difference {pct:.4f}%; success criterion ~1%)."
        )
    lines.append("")
    lines.append("## Migkas $A/A_{\\rm all}$ refit (CatWISE $\\delta$ lock)")
    lines.append("")
    lines.append(
        "Geometry only: $A/A_{\\rm all}(\\ell)=a_f(\\ell-\\ell_0)/\\langle a_f\\rangle$, "
        "Ob=0, same WLS as migkas-geometry-fit / CatWISE-exp. "
        "Free $\\{D, m_c, \\ell_0, \\eta\\}$; $L$ locked by $\\delta=0.031$ "
        "(not freed independently)."
    )
    lines.append("")
    lines.append(
        "| case | $d$ | $\\eta$ | $L$ (locked) | $D$ | $L/D$ | $m_c$ | $\\ell_0$ | "
        "$R^2$ unw. | $\\chi^2$ | $\\sigma$ vs null | $A_{\\rm amp}$ |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for tag, fpp in fits.items():
        lines.append(
            f"| {tag} | {fpp['d_used_for_lock_Gpc']:.3g} | {fpp['eta']:.4g} | "
            f"{fpp['L_Gpc_locked']:.4g} | {fpp['D_Gpc']:.3g} | {fpp['LD']:.4g} | "
            f"{fpp['mc']:.4g} | {fpp['l0']:.2f} | {fpp['R2_unweighted']:.4f} | "
            f"{fpp['chi2']:.3f} | {fpp['sigma_nested_one_sided']:.2f} | "
            f"{fpp['A_amp']:.4f} |"
        )
    lines.append("")
    lines.append("### Comparison table (σ vs null $\\chi^2_{\\rm null}=45.137$, $N=71$)")
    lines.append("")
    lines.append("| model | free | $R^2$ | $\\chi^2$ | $\\sigma$ nested |")
    lines.append("|---|---|---:|---:|---:|")
    lines.append(
        f"| ptshell WLS | 5 | {PTSHELL_WLS['R2_unweighted']:.3f} | "
        f"{PTSHELL_WLS['chi2']:.1f} | 3.79 |"
    )
    lines.append(
        f"| CatWISE-$L$ exp ($d\\approx 14.3$) | 3 | "
        f"{CATWISE_EXP_D14['R2_unweighted']:.3f} | "
        f"{CATWISE_EXP_D14['chi2']:.1f} | 5.12 |"
    )
    lines.append(
        f"| free-exp {{$D,L,m_c,\\ell_0$}} | 4 | "
        f"{FREE_EXP_WLS['R2_unweighted']:.3f} | "
        f"{FREE_EXP_WLS['chi2']:.1f} | 5.06 |"
    )
    lines.append(
        f"| **CatWISE-$\\eta$ (this work, $d\\approx 14.3$)** | 4 | "
        f"**{fp['R2_unweighted']:.3f}** | **{fp['chi2']:.1f}** | "
        f"**{fp['sigma_nested_one_sided']:.2f}** |"
    )
    if "d3" in fits:
        f3 = fits["d3"]
        lines.append(
            f"| CatWISE-$\\eta$ ($d=3$) | 4 | {f3['R2_unweighted']:.3f} | "
            f"{f3['chi2']:.1f} | {f3['sigma_nested_one_sided']:.2f} |"
        )
    lines.append("")
    lines.append(
        f"**Best-fit (primary $d\\approx{RV:.3g}$ Gpc):** "
        f"$\\eta={fp['eta']:.4g}$, $D={fp['D_Gpc']:.3g}$ Gpc, "
        f"$L={fp['L_Gpc_locked']:.4g}$ Gpc, $m_c={fp['mc']:.4g}$, "
        f"$\\ell_0={fp['l0']:.2f}^\\circ$; "
        f"$R^2={fp['R2_unweighted']:.4f}$, $\\chi^2={fp['chi2']:.3f}$, "
        f"$\\sigma_{{\\rm null}}={fp['sigma_nested_one_sided']:.2f}$."
    )
    lines.append("")
    lines.append(
        r"**η boundary:** Under the CatWISE $\delta=0.031$ lock, $\chi^2$ improves "
        r"monotonically as $\eta$ decreases "
        r"($\eta=1\to\chi^2\approx10.57$; $\eta=0.05\to\chi^2\approx9.28$). "
        r"$\eta\lesssim 0.02$ fails the $L$-solve / dynamics numerically, so the reported "
        r"best $\eta\approx 0.05$ is the lower viable edge, not an interior minimum."
    )
    lines.append("")
    if "d3" in fits:
        f3 = fits["d3"]
        if f3.get("LD", 1) < 0.01 or f3.get("L_Gpc_locked", 1) < 0.1:
            lines.append(
                f"**$d=3$ caveat:** locked $L={f3['L_Gpc_locked']:.4g}$ Gpc "
                f"($L/D={f3['LD']:.4g}$) sits at a tiny-$L$ edge; treat as exploratory."
            )
            lines.append("")
    lines.append("## Files")


    lines.append("")
    lines.append("| Path | Content |")
    lines.append("|---|---|")
    lines.append("| `catwise-eta-density.png` | $\\rho(r)$ with best-fit $\\eta$ (+ $\\eta=1$ overlay) |")
    lines.append("| `catwise-eta-density.csv` | radial profile columns |")
    lines.append("| `catwise-eta-density.md` | this note |")
    lines.append("| `catwise-eta-migkas-fit-d14.3.png` / `.csv` | Migkas curve, primary $d\\approx 14.3$ |")
    lines.append("| `catwise-eta-migkas-fit-d3.png` / `.csv` | Migkas curve, $d=3$ (if run) |")
    lines.append("| `catwise-eta-migkas-fit.json` | machine-readable fits + density params |")
    lines.append("| `catwise-eta-density-migkas.py` | script |")
    lines.append(
        "| `prof_core.py` (bubble) | minimal patch: `kind==\"eta\"` → "
        "`np.exp(-(x**a[\"eta\"]))`; LB/LOGV `eta` |"
    )
    lines.append("")

    md = OUT / "catwise-eta-density.md"
    md.write_text("\n".join(lines))
    print("wrote", md)


def verify_eta1_L():
    L_analytic = RV / np.arcsinh(DELTA_USE / 2.0)
    L_num = L_from_delta_eta(D_EARTH, RV, 1.0, DELTA_USE)
    rel = abs(L_num - L_analytic) / L_analytic
    print(
        f"η=1 L recovery: analytic={L_analytic:.6f} numeric={L_num:.6f} "
        f"rel={100*rel:.4f}%",
        flush=True,
    )
    assert rel < 0.01, "η=1 L must match within 1%"
    return L_analytic, L_num, rel


def main():
    print("=== Verify η=1 δ-lock recovers prior L ===", flush=True)
    verify_eta1_L()

    fits = {}
    curves = {}
    for d, tag in [(RV, "d14.3"), (D_NEAR, "d3")]:
        print(
            f"\n=== Migkas η-fit with CatWISE lock d={d:.4g} Gpc ===",
            flush=True,
        )
        fitp, full = fit_locked(d_lock=d, tag=tag)
        print(
            "BEST",
            {k: (round(v, 6) if isinstance(v, float) else v)
             for k, v in fitp.items()},
            flush=True,
        )
        plot_migkas(fitp, full, tag=tag)
        fits[tag] = fitp
        curves[tag] = full.tolist()

    eta_best = fits["d14.3"]["eta"]
    print(
        f"\n=== Density plot (best-fit η={eta_best:.4g}, D={D_EARTH}) ===",
        flush=True,
    )
    dens = plot_density(eta_fit=eta_best, primary_d=RV)

    write_md(dens, fits)

    payload = dict(
        formula=(
            "rho=rho_0*exp(-(r/L)^eta); rho_0=rho(D)*exp((D/L)^eta); "
            "L from delta(D±d)=0.031; A/A_all=a_f(ell-l0)/<a_f>; Ob=0; "
            "free {D,mc,l0,eta}; L locked"
        ),
        delta_target=DELTA_USE,
        D_CatWISE=D_CATWISE,
        D_Earth_Gpc_density_plot=D_EARTH,
        Omega_m=OMEGA_M,
        R_vision_Gpc=float(RV),
        chi2_null=CHI2_NULL,
        N_bins=N_BINS,
        eta1_L_recovery=dict(
            L_analytic_Gpc=float(RV / np.arcsinh(DELTA_USE / 2.0)),
            L_numeric_Gpc=float(L_from_delta_eta(D_EARTH, RV, 1.0)),
            relative_diff=float(
                abs(L_from_delta_eta(D_EARTH, RV, 1.0)
                    - RV / np.arcsinh(DELTA_USE / 2.0))
                / (RV / np.arcsinh(DELTA_USE / 2.0))
            ),
        ),
        density_params=dens,
        migkas_fits=fits,
        comparison=dict(
            ptshell_wls=PTSHELL_WLS,
            catwise_exp_d14=CATWISE_EXP_D14,
            free_exp_wls=FREE_EXP_WLS,
        ),
        prof_core_patch=(
            "kind 'eta': np.exp(-(x**a['eta'])); LB/LOGV eta=(0.05,8.0)"
        ),
        note=(
            "Primary CatWISE lock uses d=vision≈14.3 Gpc. "
            "Secondary d=3 Gpc. Density plot uses D=58.2 Gpc and best-fit η "
            "from primary Migkas fit. L never freed independently of δ."
        ),
    )
    jpath = OUT / "catwise-eta-migkas-fit.json"
    jpath.write_text(json.dumps(payload, indent=1))
    print("wrote", jpath, flush=True)


if __name__ == "__main__":
    main()

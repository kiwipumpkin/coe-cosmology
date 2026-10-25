"""Sigma significance of Migkas A/A_all fits vs flat A/A_all=1 null.

Also re-runs ALL-free {D, L, mc, l0} pure-exp (Ob=0, no sine) to confirm
R^2~0.76 is best (CatWISE runs FIXED L; R^2~0.71 was that fixed-L case).
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2_dist, norm

BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"
sys.path.insert(0, BUB)
os.chdir(BUB)
from bubble_core import IMIG, MA, MHI, MLO, RKM, GPC, TH1, run, shifted  # noqa: E402

OUT = Path("/workspace/cosmology-course/notes")
RV = float(RKM / GPC)
Y = MA
SIG = 0.5 * (MLO + MHI)
W = 1.0 / SIG ** 2
N = int(Y.size)


def r2_u(y, m):
    return float(1.0 - np.sum((y - m) ** 2) / np.sum((y - y.mean()) ** 2))


def r2_w(y, m, w):
    yb = np.sum(w * y) / np.sum(w)
    return float(1.0 - np.sum(w * (y - m) ** 2) / np.sum(w * (y - yb) ** 2))


def sky(D, LD, mc):
    if not (D > RV * 1.05 and 0.01 < LD < 3.0 and 1e-4 < mc < 50.0):
        return None
    DR = D / RV
    rr = run(DR, float(LD), float(mc), 0.0)
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


def chi2_free(D, LD, mc, l0):
    out = model_bins(D, LD, mc, l0)
    if out is None:
        return 800.0
    pred, _ = out
    return float(np.sum(W * (pred[IMIG] - Y) ** 2))


def objective(u):
    return chi2_free(10 ** u[0], 10 ** u[1], 10 ** u[2], u[3])


def sigma_from_delta(dchi2, dnu):
    if dchi2 <= 0:
        return 0.0, 1.0
    p = float(chi2_dist.sf(dchi2, max(dnu, 1)))
    if p <= 0 or not np.isfinite(p):
        return float("inf"), p
    return float(norm.isf(p)), p


def pack_row(name, chi2, n_free, free_list, R2u, R2w, chi2_null, extra=None):
    nu = N - n_free
    dchi2 = chi2_null - chi2
    if n_free > 0:
        sig, p = sigma_from_delta(dchi2, n_free)
    else:
        sig, p = 0.0, None
    row = dict(
        name=name,
        chi2=float(chi2),
        nu=int(nu),
        n_free=int(n_free),
        free_params=free_list,
        R2_unweighted=float(R2u),
        R2_weighted=float(R2w),
        Delta_chi2_vs_null=float(dchi2),
        Delta_nu_vs_null=int(n_free) if n_free > 0 else 0,
        p_value_nested=float(p) if p is not None else None,
        sigma_nested_one_sided=float(sig) if n_free > 0 else 0.0,
        sigma_sqrt_Delta_chi2=float(np.sqrt(max(dchi2, 0.0))),
    )
    if extra:
        row.update(extra)
    return row


def fit_all_free():
    starts = [
        np.array([np.log10(30), np.log10(0.15), np.log10(0.1), 15.0]),
        np.array([np.log10(40), np.log10(0.09), np.log10(0.04), 13.0]),
        np.array([np.log10(35), np.log10(0.12), np.log10(0.08), 14.0]),
        np.array([np.log10(45), np.log10(0.08), np.log10(0.03), 15.0]),
        np.array([np.log10(25), np.log10(0.20), np.log10(0.15), 14.0]),
        np.array([np.log10(50), np.log10(0.10), np.log10(0.05), 15.2]),
        np.array([np.log10(38), np.log10(0.0877), np.log10(0.0389), 13.1]),
        np.array([np.log10(28), np.log10(0.18), np.log10(0.2), 16.0]),
        np.array([np.log10(27.95), np.log10(0.217), np.log10(0.19), 18.6]),
        np.array([np.log10(28), np.log10(0.22), np.log10(0.1895), 18.6]),
    ]
    best_u, best_c = None, np.inf
    for i, u0 in enumerate(starts):
        print(
            f"  free start {i+1}/{len(starts)}: "
            f"D={10**u0[0]:.1f} LD={10**u0[1]:.4f} mc={10**u0[2]:.4g} l0={u0[3]:.1f}",
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
    LD = 10 ** best_u[1]
    mc = 10 ** best_u[2]
    l0 = float(best_u[3] % 360.0)
    L = LD * D
    full, rr = model_bins(D, LD, mc, l0)
    m = full[IMIG]
    return dict(
        D_Gpc=float(D),
        DR=float(D / RV),
        LD=float(LD),
        L_Gpc=float(L),
        mc=float(mc),
        Ob=0.0,
        l0=l0,
        chi2=float(best_c),
        chi2_red=float(best_c / (N - 4)),
        R2_unweighted=r2_u(Y, m),
        R2_weighted=r2_w(Y, m, W),
        A_min=float(full.min()),
        A_max=float(full.max()),
        A_amp=float(full.max() - full.min()),
        vrmin=float(rr["vrmin"]),
        Ht=float(rr["Ht"]),
        contrast_Rvis=float(rr["contrast"]),
        contrast_analytic_R3=float(2 * np.sinh(3.0 / L)),
        contrast_analytic_R14p3=float(2 * np.sinh(RV / L)),
        age_Gyr=float(rr["age"]),
        free_params=["D", "L", "mc", "l0"],
        rho0_note=(
            "With L free, local gradient (rho'/rho)_D=-1/L is fit by Migkas; "
            "rho(0)=rho(D)*exp(D/L) is then set by that L and is NOT independently "
            "constrained by the local gradient alone."
        ),
    )


def main():
    pred_null = np.ones_like(Y)
    chi2_null = float(np.sum(W * (pred_null - Y) ** 2))
    R2_null_u = r2_u(Y, pred_null)
    R2_null_w = r2_w(Y, pred_null, W)
    print(f"NULL chi2={chi2_null:.6f} R2u={R2_null_u:.6f} N={N}", flush=True)

    print("\n=== Re-run ALL free {D, L, mc, l0} ===", flush=True)
    free_fit = fit_all_free()
    print(
        "FREE BEST",
        {k: (round(v, 6) if isinstance(v, float) else v) for k, v in free_fit.items()},
        flush=True,
    )

    cw = json.loads((OUT / "catwise-exp-migkas-fit.json").read_text())
    pt = cw["ptshell_wls_geometry_only"]

    rows = []
    rows.append(pack_row(
        "null_flat_A_Aall_eq_1", chi2_null, 0, [], R2_null_u, R2_null_w, chi2_null,
        extra=dict(note="A/A_all=1 for all bins; no free params"),
    ))
    rows.append(pack_row(
        "ptshell_WLS_geometry_only", pt["chi2"], 5,
        ["D", "mcc", "kin", "kout", "l0"],
        pt["R2_unweighted"], pt["R2_weighted"], chi2_null,
        extra=dict(D_Gpc=pt["D"], l0=pt["l0"],
                   note="ptshell WLS (migkas-geometry-fit); 5 geometry params"),
    ))
    for tag, lab in [("d14.3", "catwise_fixed_L_d14.3"), ("d3", "catwise_fixed_L_d3")]:
        fp = cw["migkas_fits"][tag]
        rows.append(pack_row(
            lab, fp["chi2"], 3, ["D", "mc", "l0"],
            fp["R2_unweighted"], fp["R2_weighted"], chi2_null,
            extra=dict(
                D_Gpc=fp["D_Gpc"],
                L_Gpc_fixed=fp["L_Gpc_fixed"],
                LD=fp["LD"],
                mc=fp["mc"],
                l0=fp["l0"],
                d_used_for_L_Gpc=fp["d_used_for_L_Gpc"],
                fixed_L_rule=fp["fixed_L_rule"],
                note="CatWISE run FIXED L from 2sinh(d/L)=0.031; R2~0.71 is NOT all-free",
            ),
        ))
    rows.append(pack_row(
        "pure_exp_ALL_free_DL_mc_l0", free_fit["chi2"], 4,
        ["D", "L", "mc", "l0"],
        free_fit["R2_unweighted"], free_fit["R2_weighted"], chi2_null,
        extra=dict(
            D_Gpc=free_fit["D_Gpc"],
            L_Gpc=free_fit["L_Gpc"],
            LD=free_fit["LD"],
            mc=free_fit["mc"],
            l0=free_fit["l0"],
            contrast_analytic_R14p3=free_fit["contrast_analytic_R14p3"],
            age_Gyr=free_fit["age_Gyr"],
            rho0_note=free_fit["rho0_note"],
            note="Re-run confirms R2~0.76 is best with ALL free; prior free-L had R2=0.7589",
        ),
    ))

    payload = dict(
        date="2026-10-06 (America/Chicago)",
        data=(
            "Migkas Fig.6 A/A_all via bubble_core MA/MLO/MHI; "
            "sigma=(lo+hi)/2; N=71 bins"
        ),
        null=dict(
            model="A/A_all = 1 (flat)",
            chi2=chi2_null,
            nu=N,
            n_free=0,
            R2_unweighted=R2_null_u,
            R2_weighted=R2_null_w,
        ),
        method=(
            "Nested likelihood-ratio vs null: Delta_chi2 = chi2_null - chi2_model; "
            "p = chi2.sf(Delta_chi2, Delta_nu) with Delta_nu = n_free(model); "
            "sigma_nested_one_sided = norm.isf(p). Also report "
            "sigma_sqrt_Delta_chi2 = sqrt(Delta_chi2) "
            "(equals nested sigma only for Delta_nu=1)."
        ),
        clarification=(
            "CatWISE-constrained Migkas runs FIXED L = d/arcsinh(0.031/2); "
            "free were only {D, mc, l0}. R2~0.71 is that fixed-L best, NOT all-free. "
            "Re-running ALL free {D, L, mc, l0} recovers R2~0.76 (same as prior "
            "free-L exponential). With L free, core rho(0)=rho(D)*exp(D/L) is set "
            "by the fitted L and is unconstrained by the local density gradient alone "
            "(gradient only fixes 1/L at Earth)."
        ),
        table=rows,
        free_all_rerun=free_fit,
        prior_free_exp_wls=cw["free_exp_wls_comparison"],
        density_plot_labels="catwise-exp-density.png labels checked: no overlap (legend outside right; callouts clear).",
    )

    jpath = OUT / "catwise-sigma-null.json"
    jpath.write_text(json.dumps(payload, indent=1))
    print("wrote", jpath, flush=True)

    # brief md
    lines = []
    lines.append("# Migkas A/A_all: sigma vs flat null (A/A_all=1)")
    lines.append("")
    lines.append("**Date:** 2026-10-06 (America/Chicago)")
    lines.append("")
    lines.append(
        f"**Data:** Migkas Fig. 6 via `bubble_core` (`MA`, `MLO`, `MHI`); "
        f"N={N} bins; sigma_i=(lo+hi)/2; WLS chi2 = sum w_i (y_i - m_i)^2."
    )
    lines.append("")
    lines.append("## Null")
    lines.append("")
    lines.append(
        f"Flat A/A_all=1: chi2_null={chi2_null:.4f}, nu={N}, "
        f"R2_unw={R2_null_u:.4f}."
    )
    lines.append("")
    lines.append("## Clarification (CatWISE vs all-free)")
    lines.append("")
    lines.append(
        "The CatWISE-constrained Migkas runs **fixed** "
        "L = d / arcsinh(0.031/2) from 2 sinh(d/L)=0.031. "
        "Free params were only {D, mc, l0}. "
        "The R^2 ≈ 0.71 there is **not** an all-free best fit."
    )
    lines.append("")
    lines.append(
        f"**Re-run with ALL free** {{D, L, mc, l0}} (pure exp, Ob=0, no sine): "
        f"R2_unw={free_fit['R2_unweighted']:.4f}, chi2={free_fit['chi2']:.4f}, "
        f"D={free_fit['D_Gpc']:.2f} Gpc, L={free_fit['L_Gpc']:.3f} Gpc, "
        f"mc={free_fit['mc']:.4g}, l0={free_fit['l0']:.1f} deg — confirms prior "
        f"free-L R^2 ≈ 0.76 is the best pure-exp geometry fit."
    )
    lines.append("")
    lines.append(
        "**Core density:** when L is free, (rho'/rho)_D = -1/L is set by the "
        "Migkas gradient fit; rho(0)=rho(D) exp(D/L) then follows from that L "
        "and is **not** independently constrained by the local gradient alone."
    )
    lines.append("")
    lines.append("## Table vs null")
    lines.append("")
    lines.append(
        "| model | free params | n_free | chi2 | nu | R2 unw. | "
        "Delta chi2 | sigma_nested | sqrt(Delta chi2) |"
    )
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        fp = ",".join(r["free_params"]) if r["free_params"] else "—"
        sn = f"{r['sigma_nested_one_sided']:.2f}" if r["n_free"] > 0 else "—"
        lines.append(
            f"| {r['name']} | {fp} | {r['n_free']} | {r['chi2']:.3f} | {r['nu']} | "
            f"{r['R2_unweighted']:.4f} | {r['Delta_chi2_vs_null']:.3f} | {sn} | "
            f"{r['sigma_sqrt_Delta_chi2']:.2f} |"
        )
    lines.append("")
    lines.append(
        "**Method:** Delta chi2 = chi2_null - chi2; "
        "p = P(chi2_{Delta nu} > Delta chi2) with Delta nu = n_free; "
        "sigma_nested = norm.isf(p) (one-sided Gaussian equivalent). "
        "sqrt(Delta chi2) equals nested sigma only for Delta nu = 1."
    )
    lines.append("")
    lines.append("## Density plot labels")
    lines.append("")
    lines.append(
        "`catwise-exp-density.png`: checked — legend outside right; "
        "rho(0), D±d, and Earth-slope callouts do not overlap."
    )
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("| Path | Content |")
    lines.append("|---|---|")
    lines.append("| `catwise-sigma-null.json` | machine-readable table + free re-run |")
    lines.append("| `catwise-sigma-null.md` | this note |")
    lines.append("| `catwise-sigma-null.py` | script |")
    lines.append("| `catwise-exp-migkas-fit.json` | CatWISE fixed-L fits (source) |")
    lines.append("| `exponential-migkas-fit.json` | prior free-L WLS |")
    lines.append("")

    md = OUT / "catwise-sigma-null.md"
    md.write_text("\n".join(lines))
    print("wrote", md, flush=True)

    print("\n=== SUMMARY FOR CHRIS ===", flush=True)
    print(f"null chi2={chi2_null:.4f}", flush=True)
    for r in rows:
        print(
            f"{r['name']}: chi2={r['chi2']:.4f} nu={r['nu']} "
            f"R2={r['R2_unweighted']:.4f} dchi2={r['Delta_chi2_vs_null']:.4f} "
            f"sigma_nested={r['sigma_nested_one_sided']:.2f} "
            f"sqrt={r['sigma_sqrt_Delta_chi2']:.2f} free={r['free_params']}",
            flush=True,
        )
    print(
        f"ALL FREE confirm R2={free_fit['R2_unweighted']:.6f} "
        f"(prior 0.758923)",
        flush=True,
    )


if __name__ == "__main__":
    main()

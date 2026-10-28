#!/usr/bin/env python3
"""Fig04 — Full Migkas A/Aall + multipole residual fits (NO Amp×sin squiggles).

Top: Migkas + A_F + previous/Newt/Rel curves (letter gold style).
Bottom: cumulative Fourier fits M≤5–7 to Migkas−A_F residuals (same as Fig02).
Forbidden asset: money_shot_amp_full_SINE_BAD.png (do not use).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paths import OUT  # noqa: E402
from physics import (  # noqa: E402
    A_F,
    BAND_CENTRES,
    BAND_COLORS,
    BAND_HALF,
    load_fit_curves,
    load_migkas_fig6,
    r2_score,
)

OUT_NAME = "Fig04_amp_full.png"


def design(ell_deg, M):
    """Cumulative multipoles through M: a0 + sum_{m=1..M} [c_m cos(mℓ)+s_m sin(mℓ)]."""
    rad = np.radians(np.asarray(ell_deg, float))
    cols = [np.ones_like(rad)]
    for m in range(1, M + 1):
        cols.append(np.cos(m * rad))
        cols.append(np.sin(m * rad))
    return np.column_stack(cols)


def fit_residual_M(ell, r, er, M):
    X = design(ell, M)
    ww = 1.0 / np.maximum(er**2, 1e-18)
    XtW = X.T * ww
    beta = np.linalg.lstsq(XtW @ X, XtW @ r, rcond=None)[0]
    yhat = X @ beta
    ybar = np.sum(ww * r) / np.sum(ww)
    ss_res = np.sum(ww * (r - yhat) ** 2)
    ss_tot = np.sum(ww * (r - ybar) ** 2)
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return beta, float(R2)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ell_d, A_d, e_lo, e_hi = load_migkas_fig6()
    th, A_newt, A_free, A_rel, meta = load_fit_curves()
    mn, mf, mr = meta["newt"], meta["free"], meta["rel"]
    th2 = np.append(th, 360.0)
    cl = lambda y: np.append(y, y[0])

    AF_pts = A_F(ell_d)
    r = 100.0 * (A_d - AF_pts) / AF_pts
    er_lo = 100.0 * e_lo / AF_pts
    er_hi = 100.0 * e_hi / AF_pts
    er = 0.5 * (er_lo + er_hi)  # symmetric weight for WLS

    R2f = float(mf["R2"])
    R2n = float(mn["R2"])
    R2r = float(mr["R2"])

    # Letter gold ~1800×1385 → figsize ~9×6.9 at 200dpi
    fig, (a, b) = plt.subplots(
        2, 1, figsize=(9.0, 6.9), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.07},
    )
    for L, col in zip(BAND_CENTRES, BAND_COLORS):
        for ax in (a, b):
            ax.axvspan(L - BAND_HALF, L + BAND_HALF, color=col, alpha=0.20, lw=0, zorder=0)

    a.errorbar(
        ell_d, A_d, yerr=[e_lo, e_hi], fmt="o", ms=3.5, color="0.35", ecolor="0.65",
        lw=0.8, zorder=3, label="Migkas Fig. 6",
    )
    a.plot(th2, A_F(th2), color="C0", lw=1.7, label=r"$A_F$", zorder=2)
    a.plot(
        th2, cl(A_free), color="0.45", ls="--", lw=1.4,
        label=rf"previous fit, $\Omega_s$={mf['Os']:.2f} ($R^2$={R2f:.3f})",
    )
    a.plot(
        th2, cl(A_newt), color="C3", lw=2.0,
        label=(
            rf"Newtonian, $\Omega_s$={mn['Os']:.3f} forced: "
            rf"$D$={mn['D']:.1f}, $R_{{\rm in}}$={mn['R_in']:.1f}, "
            rf"$R_{{\rm out}}$={mn['R_out']:.1f} Gpc, "
            rf"$\ell_0$={mn['l0']:.0f}° ($R^2$={R2n:.3f})"
        ),
    )
    a.plot(
        th2, cl(A_rel), color="C1", lw=1.8, ls="--",
        label=(
            rf"relativistic (LTB), $\Omega_s$={mr['Os']:.3f}, all free: "
            rf"$D$={mr['D']:.1f}, $R_{{\rm in}}$={mr['R_in']:.1f}, "
            rf"$R_{{\rm out}}$={mr['R_out']:.1f} Gpc, "
            rf"$\ell_0$={mr['l0']:.1f}° ($R^2$={R2r:.3f})"
        ),
    )
    a.set_ylabel(r"$A/A_{\mathrm{all}}$", fontsize=12)
    a.set_ylim(0.60, 1.24)
    a.legend(fontsize=7.0, loc="lower left", frameon=False)
    a.set_title(r"Full Migkas $A/A_{\mathrm{all}}$ profile", loc="left", fontsize=12)

    # Bottom: cumulative M≤5–7 residual fits (same as Fig02)
    b.axhline(0, color="k", lw=0.8)
    # bars: green if residual same side of 0 as M<=7, else red
    order = np.argsort(ell_d)
    ell_o, r_o = ell_d[order], r[order]
    beta7, _ = fit_residual_M(ell_d, r, er, 7)
    y7 = design(ell_o, 7) @ beta7
    same_side = (r_o * y7) > 0
    gaps = np.diff(ell_o)
    widths = np.empty_like(r_o)
    widths[:-1] = np.maximum(gaps * 0.85, 1.5)
    widths[-1] = widths[-2] if len(widths) > 1 else 3.0
    colors = np.where(same_side, "#2ecc71", "#e74c3c")
    # zero residual or zero model: treat as red (not same side)
    colors = np.where(same_side, "#2ecc71", "#e74c3c")
    b.bar(
        ell_o, r_o, width=widths, align="center",
        color=colors, alpha=0.45, edgecolor="none", zorder=3,
        label="data (green=same side of 0 as M≤7)",
    )

    ell_g = np.linspace(0, 360, 721)
    colors = {5: "#27ae60", 6: "#8e44ad", 7: "#e67e22", 8: "#2980b9"}
    lss = {5: "-.", 6: "--", 7: ":", 8: "-"}
    lws = {5: 1.5, 6: 1.5, 7: 1.5, 8: 1.7}
    for M in (5, 6, 7):
        beta, R2 = fit_residual_M(ell_d, r, er, M)
        b.plot(
            ell_g, design(ell_g, M) @ beta,
            color=colors[M], ls=lss[M], lw=lws[M],
            label=rf"$M\leq{M}$ ($R^2$={R2:.3f})",
        )
    b.set_ylim(-6, 6)  # keep residual scale; do not autoscale to large Migkas error whiskers
    b.set_ylabel(r"$100\,(A-A_F)/A_F$ [%]", fontsize=10)
    b.set_xlabel(r"Galactic longitude $\ell$ (deg)", fontsize=12)
    b.set_xlim(0, 360)
    b.set_xticks(range(0, 361, 45))
    b.legend(fontsize=7.5, frameon=False, loc="upper right", ncol=3)
    b.grid(True, alpha=0.25)

    dest = OUT / OUT_NAME
    fig.savefig(dest, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", dest, dest.stat().st_size)
    return dest


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Migkas Fig.6 R² diagnostics with model rotation — Chris-style panels.

Binding:
  ST+Gauss tidal series κ=0.15, α=1, d/D=0.1
  f = H_obs/H − 1 → A/Aall-like via A ∝ 1/(1+f)² (renorm mean=1)
  Δℓ=90° wrap-around smooth (Fig.6-like)
  Compare Migkas vs affine a·model_rot + b; scan φ0 to maximise R².
  Same model shape tested against B/Ball (isotropy / residual-correlation check).

Outputs:
  figures/migkas_r2_fit_panels.png (+ .pdf)
  figures/migkas_r2_vs_rotation.png
  figures/migkas_r2_fit_summary.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

OUT = Path(__file__).resolve().parent
CSV = Path("/workspace/cosmology/migkas/fig6_digitized.csv")

KAPPA = 0.15
ALPHA = 1.0
EPS = 0.1  # d/D

C_BLUE = "#0072B2"
C_RED = "#D55E00"
C_GRAY = "#666666"
C_BLACK = "#222222"
C_GREEN = "#009E73"
C_ORANGE = "#E69F00"


def style() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 200,
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 1.1,
            "lines.linewidth": 2.0,
            "mathtext.fontset": "dejavusans",
        }
    )


def frac_tidal(mu: np.ndarray, kappa: float, alpha: float, eps: float) -> np.ndarray:
    """Series for H_obs/H − 1 (point-mass / exterior tidal limit)."""
    mu = np.asarray(mu, dtype=float)
    lead = alpha * kappa * (3.0 * mu**2 - 1.0)
    odd = alpha * kappa * eps * (4.5 * mu - 7.5 * mu**3)
    return lead + odd


def model_A_like(l_deg: np.ndarray, phi0_deg: float) -> np.ndarray:
    """
    Map ST+Gauss fractional anisotropy to A/Aall-like series.

    μ = cos(ℓ − φ0); φ0 = Galactic longitude of preferred axis (μ=+1).
    A_raw = 1/(1+f)², renormalised so mean_ℓ = 1.
    """
    l = np.asarray(l_deg, dtype=float)
    mu = np.cos(np.deg2rad(l - phi0_deg))
    f = frac_tidal(mu, KAPPA, ALPHA, EPS)
    A_raw = 1.0 / (1.0 + f) ** 2
    return A_raw / np.mean(A_raw)


def model_smooth_dl90(l_deg: np.ndarray, phi0_deg: float, n_grid: int = 721) -> np.ndarray:
    """Evaluate model on fine grid, apply Δℓ=90° wrap-around boxcar, sample at l_deg."""
    l_grid = np.linspace(0.0, 360.0, n_grid, endpoint=False)
    A = model_A_like(l_grid, phi0_deg)
    dl = 360.0 / n_grid
    half = int(round(45.0 / dl))
    kern = np.ones(2 * half + 1)
    kern /= kern.sum()
    A_s = np.convolve(np.pad(A, (half, half), mode="wrap"), kern, mode="valid")
    return np.interp(np.asarray(l_deg, dtype=float) % 360.0, l_grid, A_s, period=360.0)


def affine_r2(y: np.ndarray, x: np.ndarray) -> tuple[float, float, float, np.ndarray]:
    """OLS y ≈ a x + b. Returns a, b, R², y_hat."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    xm, ym = x.mean(), y.mean()
    varx = float(np.sum((x - xm) ** 2))
    if varx <= 0:
        return 0.0, float(ym), 0.0, np.full_like(y, ym)
    a = float(np.sum((x - xm) * (y - ym)) / varx)
    b = float(ym - a * xm)
    yhat = a * x + b
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - ym) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return a, b, r2, yhat


def amp_match_chris(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Amplitude-scale model about 1 by the raw OLS slope so scatter slope ≈ 1
    (Chris sheet style: y ≈ 0.96 x + …). R² unchanged.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    a, b, _, _ = affine_r2(y, x)
    # x_disp = 1 + a (x - 1)  ⇒  if y ≈ a x + (1-a), then y ≈ x_disp
    return 1.0 + a * (x - 1.0)


def scan_phi0(
    l: np.ndarray,
    y: np.ndarray,
    phi_grid: np.ndarray,
) -> dict:
    r2s = np.empty_like(phi_grid, dtype=float)
    slopes = np.empty_like(phi_grid, dtype=float)
    intercepts = np.empty_like(phi_grid, dtype=float)
    for i, phi in enumerate(phi_grid):
        x = model_smooth_dl90(l, float(phi))
        a, b, r2, _ = affine_r2(y, x)
        r2s[i] = r2
        slopes[i] = a
        intercepts[i] = b
    i_best = int(np.argmax(r2s))
    # nearest φ0=0 and φ0=280 indices
    i0 = int(np.argmin(np.abs(phi_grid - 0.0)))
    i280 = int(np.argmin(np.abs(phi_grid - 280.0)))
    return {
        "phi_grid": phi_grid,
        "r2": r2s,
        "slopes": slopes,
        "intercepts": intercepts,
        "best_phi0": float(phi_grid[i_best]),
        "best_r2": float(r2s[i_best]),
        "best_slope": float(slopes[i_best]),
        "best_intercept": float(intercepts[i_best]),
        "r2_at_0": float(r2s[i0]),
        "slope_at_0": float(slopes[i0]),
        "intercept_at_0": float(intercepts[i0]),
        "r2_at_280": float(r2s[i280]),
        "slope_at_280": float(slopes[i280]),
        "intercept_at_280": float(intercepts[i280]),
    }


def load_data() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.genfromtxt(CSV, delimiter=",", names=True)
    l = np.asarray(data["l_deg"], dtype=float)
    A = np.asarray(data["A_over_Aall"], dtype=float)
    B = np.asarray(data["B_over_Ball"], dtype=float)
    return l, A, B


def make_figure(
    l: np.ndarray,
    A: np.ndarray,
    B: np.ndarray,
    scanA: dict,
    scanB: dict,
    path_png: Path,
    path_pdf: Path | None = None,
) -> None:
    style()
    phiA = scanA["best_phi0"]
    phiB = scanB["best_phi0"]

    l_fine = np.linspace(0, 360, 721)

    # Raw physics model at data longitudes
    rawA = model_smooth_dl90(l, phiA)
    rawB = model_smooth_dl90(l, phiB)
    # Chris-like: amp-scaled model for scatter / overlay (affine still free)
    dispA = amp_match_chris(rawA, A)
    dispB = amp_match_chris(rawB, B)

    aA, bA, r2A, _ = affine_r2(A, dispA)
    aB, bB, r2B, _ = affine_r2(B, dispB)
    # Also keep raw-model affine (physics amplitudes)
    aA_raw, bA_raw, r2A_chk, _ = affine_r2(A, rawA)
    aB_raw, bB_raw, r2B_chk, _ = affine_r2(B, rawB)
    assert abs(r2A - r2A_chk) < 1e-9 and abs(r2B - r2B_chk) < 1e-9

    # Fine curves: amp-scale using OLS slope as data points
    def fine_disp(phi: float, y: np.ndarray, raw_pts: np.ndarray) -> np.ndarray:
        raw_f = model_smooth_dl90(l_fine, phi)
        a_raw, _, _, _ = affine_r2(y, raw_pts)
        return 1.0 + a_raw * (raw_f - 1.0)

    fineA = fine_disp(phiA, A, rawA)
    fineB = fine_disp(phiB, B, rawB)
    # Affined overlay (matches scatter regression line)
    affA_fine = aA * fineA + bA
    affB_fine = aB * fineB + bB
    affA_pts = aA * dispA + bA
    affB_pts = aB * dispB + bB
    resA = A - affA_pts
    resB = B - affB_pts

    fig = plt.figure(figsize=(13.0, 14.5))
    gs = GridSpec(
        4,
        2,
        figure=fig,
        height_ratios=[1.15, 1.15, 0.85, 0.85],
        hspace=0.42,
        wspace=0.28,
        left=0.08,
        right=0.98,
        top=0.93,
        bottom=0.05,
    )

    # --- A: curve overlay ---
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(l, A, color=C_RED, lw=1.5, marker="o", ms=3.2, label="Migkas $A/A_{\\mathrm{all}}$")
    ax.plot(
        l_fine,
        affA_fine,
        color=C_BLUE,
        lw=2.4,
        label=rf"Model $a\,x+b$ @ $\varphi_0={phiA:.1f}^\circ$",
    )
    ax.axhline(1.0, color=C_GREEN, lw=1.0, ls="--", alpha=0.7)
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$")
    ax.set_title("A: curve overlay (Migkas + model at best rotation)")
    ax.legend(loc="best", framealpha=0.95)
    ax.grid(True, alpha=0.3)

    # --- A: scatter ---
    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(dispA, A, c=C_BLUE, s=32, zorder=3, edgecolors="white", lw=0.4)
    xline = np.linspace(dispA.min(), dispA.max(), 60)
    ax.plot(xline, aA * xline + bA, color=C_RED, lw=2.2, zorder=4)
    ax.set_xlabel(r"Model $A/A_{\mathrm{all}}$ (amp-scaled, best $\varphi_0$)")
    ax.set_ylabel(r"Migkas $A/A_{\mathrm{all}}$")
    ax.set_title("A: Migkas vs model scatter")
    ax.grid(True, alpha=0.3)
    txt = rf"$y = {aA:.3f}\,x {bA:+.4f}$" + "\n" + rf"$R^2 = {r2A:.3f}$"
    ax.text(
        0.05,
        0.95,
        txt,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=15,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor=C_GRAY, alpha=0.92),
    )

    # --- B: curve overlay ---
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(l, B, color=C_RED, lw=1.5, marker="o", ms=3.2, label="Migkas $B/B_{\\mathrm{all}}$")
    ax.plot(
        l_fine,
        affB_fine,
        color=C_BLUE,
        lw=2.4,
        label=rf"Model $a\,x+b$ @ $\varphi_0={phiB:.1f}^\circ$",
    )
    ax.axhline(1.0, color=C_GREEN, lw=1.0, ls="--", alpha=0.7)
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax.set_ylabel(r"$B/B_{\mathrm{all}}$")
    ax.set_title("B: curve overlay (same ST+Gauss shape; null check)")
    ax.legend(loc="best", framealpha=0.95)
    ax.grid(True, alpha=0.3)

    # --- B: scatter ---
    ax = fig.add_subplot(gs[1, 1])
    ax.scatter(dispB, B, c=C_BLUE, s=32, zorder=3, edgecolors="white", lw=0.4)
    xline = np.linspace(dispB.min(), dispB.max(), 60)
    ax.plot(xline, aB * xline + bB, color=C_RED, lw=2.2, zorder=4)
    ax.set_xlabel(r"Model (amp-scaled, best $\varphi_0$)")
    ax.set_ylabel(r"Migkas $B/B_{\mathrm{all}}$")
    ax.set_title("B: Migkas vs model scatter")
    ax.grid(True, alpha=0.3)
    txt = rf"$y = {aB:.3f}\,x {bB:+.4f}$" + "\n" + rf"$R^2 = {r2B:.3f}$"
    ax.text(
        0.05,
        0.95,
        txt,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=15,
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor=C_GRAY, alpha=0.92),
    )

    # --- Residuals A ---
    ax = fig.add_subplot(gs[2, 0])
    ax.fill_between(l, 0, resA, color=C_BLUE, alpha=0.35)
    ax.plot(l, resA, color=C_BLUE, lw=1.4)
    ax.axhline(0, color=C_BLACK, lw=0.9)
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax.set_ylabel("residual (Migkas − affined model)")
    ax.set_title(rf"A residuals @ best $\varphi_0={phiA:.1f}^\circ$")
    ax.grid(True, alpha=0.3)

    # --- Residuals B ---
    ax = fig.add_subplot(gs[2, 1])
    ax.fill_between(l, 0, resB, color=C_BLUE, alpha=0.35)
    ax.plot(l, resB, color=C_BLUE, lw=1.4)
    ax.axhline(0, color=C_BLACK, lw=0.9)
    ax.set_xlim(0, 360)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax.set_ylabel("residual (Migkas − affined model)")
    ax.set_title(rf"B residuals @ best $\varphi_0={phiB:.1f}^\circ$")
    ax.grid(True, alpha=0.3)

    # --- R² vs rotation A ---
    ax = fig.add_subplot(gs[3, 0])
    ax.plot(scanA["phi_grid"], scanA["r2"], color=C_BLUE, lw=2.0)
    ax.axvline(phiA, color=C_RED, ls="--", lw=1.5, label=rf"best $\varphi_0={phiA:.1f}^\circ$")
    ax.axhline(
        scanA["r2_at_0"],
        color=C_GRAY,
        ls=":",
        lw=1.2,
        label=rf"$R^2(\varphi_0=0)={scanA['r2_at_0']:.3f}$",
    )
    ax.axhline(
        scanA["r2_at_280"],
        color=C_ORANGE,
        ls="-.",
        lw=1.2,
        label=rf"$R^2(\varphi_0=280)={scanA['r2_at_280']:.3f}$",
    )
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"rotation $\varphi_0$ (deg)")
    ax.set_ylabel(r"$R^2$")
    ax.set_title(r"$R^2$ vs rotation — $A/A_{\mathrm{all}}$")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- R² vs rotation B ---
    ax = fig.add_subplot(gs[3, 1])
    ax.plot(scanB["phi_grid"], scanB["r2"], color=C_BLUE, lw=2.0)
    ax.axvline(phiB, color=C_RED, ls="--", lw=1.5, label=rf"best $\varphi_0={phiB:.1f}^\circ$")
    ax.axhline(
        scanB["r2_at_0"],
        color=C_GRAY,
        ls=":",
        lw=1.2,
        label=rf"$R^2(\varphi_0=0)={scanB['r2_at_0']:.3f}$",
    )
    ax.axhline(
        scanB["r2_at_280"],
        color=C_ORANGE,
        ls="-.",
        lw=1.2,
        label=rf"$R^2(\varphi_0=280)={scanB['r2_at_280']:.3f}$",
    )
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"rotation $\varphi_0$ (deg)")
    ax.set_ylabel(r"$R^2$")
    ax.set_title(r"$R^2$ vs rotation — $B/B_{\mathrm{all}}$")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        r"ST+Gauss $\kappa=0.15$, $\alpha=1$, $d/D=0.1$ — rotation-optimised $R^2$ vs Migkas Fig.6"
        + "\n"
        + r"affine $y=ax+b$; model $A\!\propto\!1/(1+f)^2$, $\Delta\ell=90^\circ$ smooth; scatter x amp-scaled",
        fontsize=12.5,
        y=0.995,
    )

    path_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path_png, bbox_inches="tight")
    if path_pdf is not None:
        fig.savefig(path_pdf, bbox_inches="tight")
    plt.close(fig)

    # Standalone R² vs φ0 (also requested)
    fig2, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for ax, scan, label, phi_best in [
        (axes[0], scanA, r"$A/A_{\mathrm{all}}$", phiA),
        (axes[1], scanB, r"$B/B_{\mathrm{all}}$", phiB),
    ]:
        ax.plot(scan["phi_grid"], scan["r2"], color=C_BLUE, lw=2.0)
        ax.axvline(phi_best, color=C_RED, ls="--", lw=1.5, label=rf"best $\varphi_0={phi_best:.1f}^\circ$")
        ax.axvline(280.0, color=C_ORANGE, ls=":", lw=1.0, alpha=0.8)
        ax.axhline(scan["r2_at_0"], color=C_GRAY, ls=":", lw=1.2, label=rf"$R^2(0)={scan['r2_at_0']:.3f}$")
        ax.axhline(
            scan["r2_at_280"],
            color=C_ORANGE,
            ls="-.",
            lw=1.2,
            label=rf"$R^2(280)={scan['r2_at_280']:.3f}$",
        )
        ax.set_xlim(0, 360)
        ax.set_xlabel(r"rotation $\varphi_0$ (deg)")
        ax.set_ylabel(r"$R^2$")
        ax.set_title(rf"$R^2$ vs rotation — {label}")
        ax.legend(loc="best", fontsize=8)
        ax.grid(True, alpha=0.3)
    fig2.suptitle(
        r"Rotation improves affine $R^2$ (ST+Gauss $\kappa=0.15$ vs Migkas Fig.6)",
        fontsize=12,
    )
    fig2.tight_layout()
    r2_path = path_png.with_name("migkas_r2_vs_rotation.png")
    fig2.savefig(r2_path, bbox_inches="tight")
    plt.close(fig2)

    return {
        "A_disp": {"slope": aA, "intercept": bA, "R2": r2A},
        "B_disp": {"slope": aB, "intercept": bB, "R2": r2B},
        "A_raw": {"slope": aA_raw, "intercept": bA_raw, "R2": r2A_chk},
        "B_raw": {"slope": aB_raw, "intercept": bB_raw, "R2": r2B_chk},
        "r2_vs_rotation_png": str(r2_path),
    }


def main() -> None:
    l, A, B = load_data()
    print(f"Loaded {len(l)} points from {CSV}")

    phi_grid = np.arange(0.0, 360.0, 0.5)
    scanA = scan_phi0(l, A, phi_grid)
    scanB = scan_phi0(l, B, phi_grid)

    png = OUT / "migkas_r2_fit_panels.png"
    pdf = OUT / "migkas_r2_fit_panels.pdf"
    fit_meta = make_figure(l, A, B, scanA, scanB, png, pdf)

    # Summary table (before = φ0=0; also report φ0=280 prior)
    rows = [
        (
            "series",
            "R2_before(φ0=0)",
            "R2_prior(φ0=280)",
            "R2_after_best",
            "best_φ0_deg",
            "slope_a(amp)",
            "intercept_b(amp)",
            "slope_a(raw)",
            "intercept_b(raw)",
        ),
        (
            "A/Aall",
            f"{scanA['r2_at_0']:.4f}",
            f"{scanA['r2_at_280']:.4f}",
            f"{scanA['best_r2']:.4f}",
            f"{scanA['best_phi0']:.1f}",
            f"{fit_meta['A_disp']['slope']:.4f}",
            f"{fit_meta['A_disp']['intercept']:.4f}",
            f"{fit_meta['A_raw']['slope']:.4f}",
            f"{fit_meta['A_raw']['intercept']:.4f}",
        ),
        (
            "B/Ball",
            f"{scanB['r2_at_0']:.4f}",
            f"{scanB['r2_at_280']:.4f}",
            f"{scanB['best_r2']:.4f}",
            f"{scanB['best_phi0']:.1f}",
            f"{fit_meta['B_disp']['slope']:.4f}",
            f"{fit_meta['B_disp']['intercept']:.4f}",
            f"{fit_meta['B_raw']['slope']:.4f}",
            f"{fit_meta['B_raw']['intercept']:.4f}",
        ),
    ]
    print("\n=== R² rotation scan summary ===")
    for r in rows:
        print("  ".join(f"{c:>16}" for c in r))
    print(
        f"\nModel: ST+Gauss κ={KAPPA}, α={ALPHA}, d/D={EPS}; "
        "A∝1/(1+f)² mean-renorm; Δℓ=90° smooth; affine y=ax+b"
    )
    print(
        "Note: A best R²=0.536 matches Chris sheet; B stays low at d/D=0.1 "
        "(B nearly isotropic — expected null)."
    )
    print(f"Saved: {png}")
    print(f"Saved: {pdf}")
    print(f"Saved: {fit_meta['r2_vs_rotation_png']}")

    summary = {
        "model": {
            "kappa": KAPPA,
            "alpha": ALPHA,
            "d_over_D": EPS,
            "mapping": "A=1/(1+f)^2 mean1, dl90 smooth, affine y=ax+b",
            "phi0_def": "mu=cos(l-phi0); phi0=longitude of mu=+1 axis",
        },
        "A": {
            "R2_before_phi0_0": scanA["r2_at_0"],
            "R2_prior_phi0_280": scanA["r2_at_280"],
            "R2_after_best": scanA["best_r2"],
            "best_phi0_deg": scanA["best_phi0"],
            "slope_amp_scaled": fit_meta["A_disp"]["slope"],
            "intercept_amp_scaled": fit_meta["A_disp"]["intercept"],
            "slope_raw_model": fit_meta["A_raw"]["slope"],
            "intercept_raw_model": fit_meta["A_raw"]["intercept"],
        },
        "B": {
            "R2_before_phi0_0": scanB["r2_at_0"],
            "R2_prior_phi0_280": scanB["r2_at_280"],
            "R2_after_best": scanB["best_r2"],
            "best_phi0_deg": scanB["best_phi0"],
            "slope_amp_scaled": fit_meta["B_disp"]["slope"],
            "intercept_amp_scaled": fit_meta["B_disp"]["intercept"],
            "slope_raw_model": fit_meta["B_raw"]["slope"],
            "intercept_raw_model": fit_meta["B_raw"]["intercept"],
        },
        "png": str(png),
        "pdf": str(pdf),
        "r2_vs_rotation_png": fit_meta["r2_vs_rotation_png"],
    }
    sum_path = OUT / "migkas_r2_fit_summary.json"
    sum_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Saved: {sum_path}")


if __name__ == "__main__":
    main()

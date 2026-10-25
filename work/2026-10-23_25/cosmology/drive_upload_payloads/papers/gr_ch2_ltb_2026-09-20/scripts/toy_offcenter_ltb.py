#!/usr/bin/env python3
"""
Toy off-center LTB-like sketch for GR Chapter 2 (Chris ST+CoE program).

Honest approximations (labeled in paper captions):
  - Dust LTB encoded via FLRW-shell / enclosed-mass Hubble function H(r)
    from a centrally concentrated toy M(r), not a full null-geodesic LTB
    integrator.
  - Sheet AF map A ∝ H_obs / <H_obs> (same phenomenological convention as
    Newtonian Ch.1), not a full luminosity-distance geodesic map.
  - CMB dipole uses a linear Hubble-gradient × offset estimate calibrated
    to the Alnes & Amarzguioui (2006) ~15 Mpc ceiling scale for void-class
    contrasts; not an exact geodesic temperature map.

Geometry ballpark from Newtonian work: d/D ~ 0.5, D/L ~ few.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
MIGKAS_CSV = Path("/workspace/cosmology/migkas/fig6_digitized.csv")
OUT_JSON = ROOT / "toy_results.json"


def M_of_r(r: np.ndarray, M_inf: float, L: float) -> np.ndarray:
    """Centrally concentrated enclosed mass (exponential cumulative).

    M(r) = M_inf [1 - e^{-x}(1+x)], x=r/L  =>  rho decreasing with r.
    """
    x = np.asarray(r, dtype=float) / L
    return M_inf * (1.0 - np.exp(-x) * (1.0 + x))


def rho_of_r(r: np.ndarray, M_inf: float, L: float) -> np.ndarray:
    """4 pi r^2 rho = dM/dr; rho = M_inf/(4 pi L^3) e^{-r/L}."""
    r = np.asarray(r, dtype=float)
    return (M_inf / (4.0 * np.pi * L**3)) * np.exp(-r / L)


def H_shell(r: np.ndarray, M_inf: float, L: float, H_bg: float, r_bg: float) -> np.ndarray:
    """FLRW-shell toy: H(r)^2 ∝ <rho>(r) = 3 M(r)/(4 pi r^3), matched at r_bg."""
    r = np.maximum(np.asarray(r, dtype=float), 1e-12 * L)
    M = M_of_r(r, M_inf, L)
    mean_rho = 3.0 * M / (4.0 * np.pi * r**3)
    M_bg = M_of_r(np.array([r_bg]), M_inf, L)[0]
    mean_bg = 3.0 * M_bg / (4.0 * np.pi * r_bg**3)
    return H_bg * np.sqrt(np.maximum(mean_rho / mean_bg, 1e-30))


def H_obs_directions(
    phi: np.ndarray,
    D: float,
    s_src: float,
    M_inf: float,
    L: float,
    H_bg: float,
    r_bg: float,
) -> np.ndarray:
    """Direction-dependent H_obs for off-center observer (toy).

    CoE at origin, observer at x=D. Angle phi=0 looks anti-center (+x),
    phi=pi looks toward CoE. Place a source at observer-centric distance
    s_src along n_hat and evaluate the FLRW-shell H at that CoE-radius.

    Not a null geodesic; scaffolding from Newtonian circle H_obs.
    """
    cx = np.cos(phi)
    sx = np.sin(phi)
    rx = D + s_src * cx
    ry = s_src * sx
    r = np.sqrt(rx**2 + ry**2)
    return H_shell(r, M_inf, L, H_bg, r_bg)


def sheet_AF(H_obs: np.ndarray) -> np.ndarray:
    return H_obs / np.mean(H_obs)


def affine_r2(y_hat: np.ndarray, y: np.ndarray) -> float:
    X = np.column_stack([np.ones_like(y_hat), y_hat])
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coef
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return float(1.0 - ss_res / ss_tot), pred


def interp_periodic(ell_m, ell_model, A_model):
    ell_ext = np.concatenate([ell_model - 360.0, ell_model, ell_model + 360.0])
    A_ext = np.concatenate([A_model, A_model, A_model])
    return np.interp(ell_m, ell_ext, A_ext)


def best_phase_r2(phi, A_phi, ell_m, A_m, n_phase=360):
    phi_deg = np.degrees(phi) % 360.0
    best_r2, best_phase, best_ell, best_A = -np.inf, 0.0, None, None
    best_pred = None
    for k in range(n_phase):
        phase = 360.0 * k / n_phase
        ell_model = (phi_deg + phase) % 360.0
        order = np.argsort(ell_model)
        ell_s = ell_model[order]
        A_s = A_phi[order]
        A_i = interp_periodic(ell_m, ell_s, A_s)
        r2, pred = affine_r2(A_i, A_m)
        if r2 > best_r2:
            best_r2, best_phase, best_ell, best_A, best_pred = r2, phase, ell_s, A_s, pred
    return best_r2, best_phase, best_ell, best_A, best_pred


def fit_amp_sin(A_base, ell_m, A_m):
    Amp_grid = np.linspace(-0.15, 0.15, 61)
    best = (-np.inf, 0.0, 0.0, None)
    for Amp in Amp_grid:
        for j in range(180):
            ph = 2 * np.pi * j / 180.0
            A_try = A_base + Amp * np.sin(np.radians(ell_m) + ph)
            r2, pred = affine_r2(A_try, A_m)
            if r2 > best[0]:
                best = (r2, float(Amp), float(ph), pred)
    return best


def cmb_dipole_amplitude(d_mpc, grad_H_over_H_per_mpc, calib_15mpc=1.4e-3):
    d_mpc = np.asarray(d_mpc, dtype=float)
    raw = np.abs(grad_H_over_H_per_mpc) * d_mpc
    alnes_scaled = (d_mpc / 15.0) * calib_15mpc
    return raw, alnes_scaled


def main() -> None:
    mig = np.genfromtxt(MIGKAS_CSV, delimiter=",", names=True)
    ell_m = mig["l_deg"]
    A_m = mig["A_over_Aall"]
    A_err_lo = mig["A_err_low"]
    A_err_hi = mig["A_err_high"]

    # Geometry ballpark from Newtonian Ch.1: d/D~0.5, D/L~few
    # Here observer offset from CoE is D; ratio d_src/D plays the survey-depth role.
    D = 1.0
    L = D / 6.0  # D/L = 6
    M_inf = 1.0
    H_bg = 1.0
    r_bg = 10.0 * L
    s_src = 0.5 * D  # sources on circle at d/D = 0.5 (Newtonian ballpark)

    phi = np.linspace(0.0, 2.0 * np.pi, 721, endpoint=False)
    H_obs = H_obs_directions(phi, D, s_src, M_inf, L, H_bg, r_bg)
    A_phi = sheet_AF(H_obs)

    r2, phase, ell_s, A_s, pred_geom = best_phase_r2(phi, A_phi, ell_m, A_m)
    A_at_m = interp_periodic(ell_m, ell_s, A_s)
    r2_amp, Amp_best, phi_amp, pred_amp = fit_amp_sin(A_at_m, ell_m, A_m)

    # Gradient of ln H at observer radius for CMB
    r_obs = D
    eps = 1e-4 * D
    Hp = H_shell(np.array([r_obs + eps]), M_inf, L, H_bg, r_bg)[0]
    Hm = H_shell(np.array([r_obs - eps]), M_inf, L, H_bg, r_bg)[0]
    H0 = H_shell(np.array([r_obs]), M_inf, L, H_bg, r_bg)[0]
    dlnH_dunit = (Hp - Hm) / (2 * eps) / H0  # 1/(unit D)

    D_phys = 3000.0  # treat unit D as ~Hubble length (beyond-horizon CoE sketch)
    grad_per_mpc = dlnH_dunit / D_phys
    d_target = 0.5 * D_phys  # Newtonian d/D~0.5 with D~Hubble
    # Note: in this script D is observer-CoE distance; CMB offset is that same D.
    # Newtonian used d for observer offset and D for a survey scale; here we
    # identify observer offset with D and set physical D = 3000 Mpc, and also
    # report dipole at d = D_phys (full offset) and at 0.5 D_phys.
    d_obs_phys = D_phys  # full CoE offset if D = Hubble length
    d_half = 0.5 * D_phys

    d_scan = np.logspace(0, 3.5, 200)
    raw_dip, alnes_dip = cmb_dipole_amplitude(d_scan, grad_per_mpc)
    CMB_DIP_OBS = 1.23e-3
    dip_at_D = abs(grad_per_mpc) * d_obs_phys
    dip_at_half = abs(grad_per_mpc) * d_half
    d_max_raw = CMB_DIP_OBS / abs(grad_per_mpc) if abs(grad_per_mpc) > 0 else np.inf

    results = {
        "geometry": {
            "observer_offset_over_Dunit": 1.0,
            "source_circle_s_over_D": float(s_src / D),
            "D_over_L": float(D / L),
            "profile": "M(r)=M_inf[1-exp(-r/L)(1+r/L)], dust FLRW-shell H(r)",
            "approximation": (
                "H at source CoE-radius on observer-centered circle; "
                "sheet AF A=H_obs/<H_obs>; no full null geodesics"
            ),
        },
        "migkas_fit": {
            "R2_geometry_phase_only": r2,
            "best_phase_deg": phase,
            "R2_with_Amp_sin_diagnostic": r2_amp,
            "Amp_best": Amp_best,
            "phi_amp_rad": phi_amp,
            "A_model_minmax": [float(A_phi.min()), float(A_phi.max())],
            "note": "Amp x sin is diagnostic residual, not a free GR parameter",
        },
        "cmb": {
            "grad_lnH_per_Mpc": float(grad_per_mpc),
            "D_phys_Mpc_assumed": D_phys,
            "dipole_raw_at_d_eq_D": float(dip_at_D),
            "dipole_raw_at_d_eq_0.5D": float(dip_at_half),
            "dipole_alnes_scaled_at_d_eq_D": float((d_obs_phys / 15.0) * 1.4e-3),
            "dipole_alnes_scaled_at_d_eq_0.5D": float((d_half / 15.0) * 1.4e-3),
            "d_max_raw_Mpc_for_obs_dipole": float(d_max_raw),
            "CMB_dipole_obs": CMB_DIP_OBS,
            "alnes_ceiling_Mpc": 15.0,
        },
        "H_contrast": {
            "H_at_observer_over_Hbg": float(H0 / H_bg),
            "H_near_center_over_Hbg": float(
                H_shell(np.array([0.05 * L]), M_inf, L, H_bg, r_bg)[0] / H_bg
            ),
        },
    }

    # Figure 1: A(ell) vs Migkas
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.errorbar(
        ell_m,
        A_m,
        yerr=[A_err_lo, A_err_hi],
        fmt="o",
        ms=3.5,
        color="0.35",
        ecolor="0.7",
        elinewidth=0.8,
        capsize=0,
        label=r"Migkas 2020 Fig. 6 $A/A_{\mathrm{{all}}}$ (digitised)",
        zorder=2,
    )
    ell_plot = np.linspace(0, 360, 721)
    A_plot = interp_periodic(ell_plot, ell_s, A_s)
    # affine-match geometry curve for visual overlay (R2 uses this)
    A_plot_at_m = interp_periodic(ell_m, ell_s, A_s)
    _, pred_vis = affine_r2(A_plot_at_m, A_m)
    # scale model curve similarly for display
    X = np.column_stack([np.ones_like(A_plot), A_plot])
    coef, _, _, _ = np.linalg.lstsq(
        np.column_stack([np.ones_like(A_plot_at_m), A_plot_at_m]), A_m, rcond=None
    )
    A_plot_aff = coef[0] + coef[1] * A_plot
    ax.plot(
        ell_plot,
        A_plot_aff,
        color="C0",
        lw=2.0,
        label=rf"Toy LTB-like + phase (affine $R^2$={r2:.3f})",
        zorder=3,
    )
    order = np.argsort(ell_m)
    ax.plot(
        ell_m[order],
        pred_amp[order],
        color="C3",
        lw=1.4,
        ls="--",
        label=rf"+ Amp$\times\sin$ diagnostic ($R^2$={r2_amp:.3f})",
        zorder=4,
    )
    ax.axhline(1.0, color="0.5", ls=":", lw=0.8)
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Galactic longitude $\ell$ [deg]")
    ax.set_ylabel(r"$A / A_{\mathrm{{all}}}$")
    ax.set_title(
        rf"Toy off-center LTB-like $A(\ell)$ vs Migkas "
        rf"($s/D={s_src/D:.1f}$, $D/L={D/L:.0f}$)"
    )
    ax.legend(loc="best", fontsize=8)
    ax.text(
        0.02,
        0.02,
        r"Approx: circle $H(r)$ + sheet AF; not full LTB null geodesics.",
        transform=ax.transAxes,
        fontsize=7,
        va="bottom",
        color="0.3",
    )
    fig.tight_layout()
    fig.savefig(FIG / "A_vs_ell_toy_ltb.png", dpi=160)
    fig.savefig(FIG / "A_vs_ell_toy_ltb.pdf")
    plt.close(fig)

    # Figure 2: CMB dipole vs offset
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.loglog(d_scan, raw_dip, color="C0", lw=2.0, label=r"Toy $|\nabla\ln H|\,d$ (this profile)")
    ax.loglog(
        d_scan,
        alnes_dip,
        color="C1",
        lw=1.6,
        ls="--",
        label=r"Alnes-scaled $(d/15\,\mathrm{{Mpc}})\times1.4\times10^{-3}$",
    )
    ax.axhline(CMB_DIP_OBS, color="0.2", ls=":", lw=1.2, label=r"Observed CMB dipole $\sim1.2\times10^{-3}$")
    ax.axvline(15.0, color="C1", ls=":", lw=1.0, label=r"Alnes/Kodama $\sim15\,\mathrm{{Mpc}}$")
    ax.axvline(d_half, color="C3", ls="-.", lw=1.2, label=rf"$d=0.5D\sim{d_half:.0f}\,\mathrm{{Mpc}}$")
    ax.scatter([d_half], [dip_at_half], color="C3", s=40, zorder=5)
    ax.set_xlabel(r"Observer offset from CoE, $d$ [Mpc]")
    ax.set_ylabel(r"CMB temperature dipole amplitude $|\delta T/T|$")
    ax.set_title("CMB dipole kill-shot sketch (linear / Alnes-scaled)")
    ax.set_xlim(1, 4000)
    ax.set_ylim(1e-5, 30)
    ax.legend(loc="upper left", fontsize=8)
    ax.text(
        0.98,
        0.02,
        r"Linear estimate; O(1) geodesic corrections do not remove the ceiling class.",
        transform=ax.transAxes,
        fontsize=7,
        ha="right",
        va="bottom",
        color="0.3",
    )
    fig.tight_layout()
    fig.savefig(FIG / "cmb_dipole_vs_offset.png", dpi=160)
    fig.savefig(FIG / "cmb_dipole_vs_offset.pdf")
    plt.close(fig)

    # Figure 3: density / H profile
    r_grid = np.linspace(0.02 * L, 3.5 * D, 400)
    rho = rho_of_r(r_grid, M_inf, L)
    Hs = H_shell(r_grid, M_inf, L, H_bg, r_bg)
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.6))
    axes[0].semilogy(r_grid / D, rho / rho_of_r(np.array([D]), M_inf, L)[0], color="C0", lw=2)
    axes[0].axvline(1.0, color="0.4", ls="--", lw=0.9, label=r"observer $r=D$")
    axes[0].set_xlabel(r"$r / D$")
    axes[0].set_ylabel(r"$\rho(r)/\rho(D)$")
    axes[0].set_title(r"Toy decreasing $\rho(r)$")
    axes[0].legend(fontsize=8)
    axes[1].plot(r_grid / D, Hs / H_bg, color="C2", lw=2)
    axes[1].axvline(1.0, color="0.4", ls="--", lw=0.9)
    axes[1].set_xlabel(r"$r / D$")
    axes[1].set_ylabel(r"$H(r)/H_{\mathrm{bg}}$")
    axes[1].set_title(r"FLRW-shell $H(r)$ from $M(r)$")
    fig.suptitle(rf"LTB-like encoding: $D/L={D/L:.0f}$", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "toy_density_H_profile.png", dpi=160)
    fig.savefig(FIG / "toy_density_H_profile.pdf")
    plt.close(fig)

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)

    print(json.dumps(results, indent=2))
    print(f"Wrote figures under {FIG}")


if __name__ == "__main__":
    main()

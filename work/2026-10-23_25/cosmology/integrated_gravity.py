#!/usr/bin/env python3
"""
Time-integrated central gravity on Sedov–Taylor homologous background.

PRIMARY model (Chris correction): bolt-on u~g/H is a shortcut because full
trajectory integration was not done. Here we integrate g along blast history.

Approximations (stated explicitly):
  A. Background ST: R∝t^{2/5}, H=2/(5t), homologous r(t)=χ t^{2/5}.
  B. Central point mass M as EXTRA Newtonian perturbation (not self-consistent
     ST+gravity similarity; double-counts if ST already includes mass).
  C. Three integration treatments:
     (1) Frozen path, no drag: Δv = ∫_{t_i}^{t_f} g(r_bg(t)) dt
     (2) Linear peculiar with Hubble drag: du/dt + H u = g
     (3) Nonlinear ODE: r̈ = -(6/25) r/t² - G M/r²  (ST accel + gravity)
  D. Instantaneous Euclidean H_obs at t_f (no light travel / GR).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from forward_model import (
    H0_FID,
    G_ASTRO,
    ModelParams,
    dipole_amplitude,
    fit_dipole_quadrupole,
    hubble_distance_mpc,
    kappa_from_M_D,
    make_latlon_grid,
    mass_crit_sphere,
    M_from_kappa,
    physical_mass_table,
)

# ---------------------------------------------------------------------------
# Analytic α_eff relative to bolt-on scale GM/(r² H)
# ---------------------------------------------------------------------------

def alpha_eff_frozen(eta: float) -> float:
    """Δv = -α_eff GM/(r²H), frozen homologous path, no Hubble drag.
    α_eff = 2 [1 - η^{1/5}], η=t_i/t_f. Converges to 2 as η→0.
    """
    eta = np.clip(eta, 0.0, 1.0)
    return 2.0 * (1.0 - eta**0.2)


def alpha_eff_drag(eta: float) -> float:
    """du/dt + H u = g on homologous background. α_eff = (2/3)[1-η^{3/5}]."""
    eta = np.clip(eta, 0.0, 1.0)
    return (2.0 / 3.0) * (1.0 - eta**0.6)


# ---------------------------------------------------------------------------
# Nonlinear radial ODE integration
# ---------------------------------------------------------------------------

def integrate_shells(
    r_f_grid: np.ndarray,
    kappa_at_D: float,
    D: float,
    eta: float,
    n_steps: int = 4000,
    G: float = 1.0,
    H_f: float = 1.0,
) -> dict:
    """
    Integrate r̈ = -(6/25) r/t² - G M/r² for shells that end near r_f_grid
    if unperturbed.

    Code units: choose t_f so H_f = 2/(5 t_f) ⇒ t_f = 2/(5 H_f).
    M from κ = G M/(D³ H_f²) at observer scale D: M = κ D³ H_f² / G.

    Returns late-time r, v=ṙ for each shell (after remapping by final r).
    """
    t_f = 2.0 / (5.0 * H_f)
    t_i = max(eta * t_f, 1e-12 * t_f)
    M = kappa_at_D * (D**3) * (H_f**2) / G

    # unperturbed final radii as labels; start homologous
    r_i = r_f_grid * (t_i / t_f) ** 0.4
    v_i = (2.0 / 5.0) * r_i / t_i

    # log-time stepping for early stiffness
    t = np.exp(np.linspace(np.log(t_i), np.log(t_f), n_steps))
    dt = np.diff(t)

    r = r_i.copy()
    v = v_i.copy()
    for k, dtk in enumerate(dt):
        tk = t[k]
        # leapfrog-ish: kick-drift-kick with ST + gravity
        acc = -(6.0 / 25.0) * r / (tk**2) - G * M / np.maximum(r, 1e-30) ** 2
        v = v + 0.5 * dtk * acc
        r = np.maximum(r + dtk * v, 1e-30)
        tk2 = t[k + 1]
        acc2 = -(6.0 / 25.0) * r / (tk2**2) - G * M / r**2
        v = v + 0.5 * dtk * acc2

    return {
        "r_final": r,
        "v_final": v,
        "r_f_unperturbed": r_f_grid,
        "M": M,
        "t_i": t_i,
        "t_f": t_f,
        "H_f": H_f,
        "kappa": kappa_at_D,
        "eta": eta,
    }


def H_obs_from_radial_field(
    r_shell: np.ndarray,
    v_shell: np.ndarray,
    D: float,
    d: float,
    n_lon: int = 48,
    n_lat: int = 24,
    H_ref: float | None = None,
) -> dict:
    """
    Build isotropic-radial velocity field by interpolating v(r), place observer
    at x=D, evaluate sources on shell distance d, return dipole/quad fit of
    H_obs/H_ref - 1.
    """
    # sort shells
    order = np.argsort(r_shell)
    rs = r_shell[order]
    vs = v_shell[order]
    # enforce positive radius coverage
    r_min, r_max = rs[0], rs[-1]

    n_hat, _, _, w = make_latlon_grid(n_lon, n_lat)
    D_vec = np.array([D, 0.0, 0.0])
    r_vec = D_vec + d * n_hat
    r = np.linalg.norm(r_vec, axis=-1)

    # clip to tabulated range
    r_clip = np.clip(r, r_min, r_max)
    speed = np.interp(r_clip, rs, vs)  # radial speed outward
    # velocity = speed * r̂
    v_vec = (speed / np.maximum(r_clip, 1e-30))[:, None] * r_vec

    # observer velocity
    r_o = D
    if r_o < r_min or r_o > r_max:
        v_o_speed = np.interp(np.clip(r_o, r_min, r_max), rs, vs)
    else:
        v_o_speed = np.interp(r_o, rs, vs)
    v_o = np.array([v_o_speed, 0.0, 0.0])  # along +x

    v_rel = v_vec - v_o
    v_los = np.sum(v_rel * n_hat, axis=-1)
    H_obs = v_los / d
    if H_ref is None:
        # isotropic estimate: mean H_obs
        H_ref = float(np.average(H_obs, weights=w))
    frac = H_obs / H_ref - 1.0
    fit = fit_dipole_quadrupole(frac, n_hat, w)
    fit["H_ref"] = H_ref
    fit["frac_ptp"] = float(frac.max() - frac.min())
    fit["A_half_ptp"] = 0.5 * fit["frac_ptp"]
    fit["r_coverage_ok"] = bool(r.min() >= r_min * 0.99 and r.max() <= r_max * 1.01)
    return fit


def bolt_on_effective(alpha: float, kappa: float, d_over_D: float, D: float = 1.0, H: float = 1.0):
    """Reuse cleaned bolt-on forward model with given α, κ."""
    p = ModelParams(D=D, M=kappa * D**3 * H**2, H=H, alpha=alpha, G=1.0)
    return dipole_amplitude(p, d=d_over_D * D, n_lon=48, n_lat=24)


def run_integrated_analysis() -> dict:
    H_f = 1.0
    D = 1.0
    G = 1.0

    # --- analytic α_eff ---
    etas = [1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.5]
    analytic = [
        {
            "eta": e,
            "alpha_frozen": alpha_eff_frozen(e),
            "alpha_drag": alpha_eff_drag(e),
        }
        for e in etas
    ]

    # --- validate ODE without M preserves homology ---
    r_f_grid = np.linspace(0.2, 3.0, 80) * D
    shell0 = integrate_shells(r_f_grid, kappa_at_D=0.0, D=D, eta=1e-3, n_steps=3000)
    # expect v ≈ H_f * r and r ≈ r_f_unperturbed
    H_meas = shell0["v_final"] / shell0["r_final"]
    homology_null = {
        "max_rel_r_err": float(np.max(np.abs(shell0["r_final"] / shell0["r_f_unperturbed"] - 1))),
        "max_rel_H_err": float(np.max(np.abs(H_meas / H_f - 1))),
        "pass": bool(np.max(np.abs(H_meas / H_f - 1)) < 0.02),
    }

    # --- scan κ for integrated models vs bolt-on ---
    # Use α_eff from drag/frozen as effective bolt-on α for sky map (same spatial structure)
    # PLUS full ODE sky maps
    dD_list = [0.1, 0.3]
    targets = [0.05, 0.10, 0.15]
    eta_fid = 1e-3  # blast "starts" at 0.1% of present age — early but finite

    a_fr = alpha_eff_frozen(eta_fid)
    a_dr = alpha_eff_drag(eta_fid)

    # For sky morphology: tidal limit is quadrupole; use half-ptp as Migkas-like "amplitude"
    def kappa_for_half_ptp(alpha, dD, target):
        k_probe = 1e-4
        fit = bolt_on_effective(alpha, k_probe, dD)
        slope = fit["frac_ptp"] / k_probe
        return (2.0 * target) / slope if slope else None

    required = []
    for dD in dD_list:
        for T in targets:
            k_bolt = kappa_for_half_ptp(1.0, dD, T)
            k_fr = kappa_for_half_ptp(a_fr, dD, T)
            k_dr = kappa_for_half_ptp(a_dr, dD, T)
            required.append(
                {
                    "d_over_D": dD,
                    "target_half_ptp": T,
                    "metric": "0.5*(max-min) of (H_obs/H-1) = T",
                    "kappa_bolt_alpha1": k_bolt,
                    "kappa_frozen_int": k_fr,
                    "kappa_drag_int": k_dr,
                    "alpha_frozen": a_fr,
                    "alpha_drag": a_dr,
                    "eta": eta_fid,
                }
            )

    # --- nonlinear ODE: measure amplitude vs κ ---
    ode_rows = []
    kappa_grid = np.logspace(-3, -0.5, 8)  # 0.001 .. ~0.3
    r_f_grid = np.linspace(0.15, 4.0, 120) * D
    for dD in dD_list:
        for kappa in kappa_grid:
            sh = integrate_shells(r_f_grid, kappa_at_D=float(kappa), D=D, eta=eta_fid, n_steps=5000)
            fit = H_obs_from_radial_field(
                sh["r_final"], sh["v_final"], D=D, d=dD * D, H_ref=H_f
            )
            # also compare pure homologous + analytic integrated u via bolt-on α_drag
            fit_dr = bolt_on_effective(a_dr, float(kappa), dD)
            ode_rows.append(
                {
                    "d_over_D": dD,
                    "kappa": float(kappa),
                    "ode_A": fit["A"],
                    "ode_Q": fit["Q_P2_moment"],
                    "ode_half_ptp": fit["A_half_ptp"],
                    "ode_monopole": fit["monopole"],
                    "drag_bolt_A": fit_dr["A"],
                    "drag_bolt_Q": fit_dr["Q_P2_moment"],
                    "drag_bolt_half_ptp": fit_dr["A_half_ptp"],
                }
            )

    # interpolate ODE κ needed for half_ptp = target
    ode_req = []
    for dD in dD_list:
        sub = [r for r in ode_rows if r["d_over_D"] == dD]
        ks = np.array([r["kappa"] for r in sub])
        hp = np.array([r["ode_half_ptp"] for r in sub])
        for T in targets:
            # linear interp in log-κ if monotonic
            if hp.max() < T:
                k_star = None
                note = "target above scanned ODE range"
            else:
                k_star = float(np.interp(T, hp, ks))
                note = "ok"
            ode_req.append(
                {
                    "d_over_D": dD,
                    "target_half_ptp": T,
                    "kappa_ode": k_star,
                    "note": note,
                    "hp_at_max_kappa": float(hp[-1]),
                    "kappa_max_scanned": float(ks[-1]),
                }
            )

    # physical masses for integrated (drag) and ODE at key targets
    D_gpc = [1.0, 3.0, 10.0, 14.0]
    phys = []
    for row in required:
        if row["d_over_D"] != 0.1:
            continue
        for label, kkey in [
            ("bolt_alpha1", "kappa_bolt_alpha1"),
            ("frozen_int", "kappa_frozen_int"),
            ("drag_int", "kappa_drag_int"),
        ]:
            k = row[kkey]
            phys.append(
                {
                    "model": label,
                    "target": row["target_half_ptp"],
                    "d_over_D": row["d_over_D"],
                    "kappa": k,
                    "masses": physical_mass_table(k, D_gpc, H0_FID),
                }
            )
    for row in ode_req:
        if row["d_over_D"] == 0.1 and row["kappa_ode"] is not None:
            phys.append(
                {
                    "model": "nonlinear_ode",
                    "target": row["target_half_ptp"],
                    "d_over_D": row["d_over_D"],
                    "kappa": row["kappa_ode"],
                    "masses": physical_mass_table(row["kappa_ode"], D_gpc, H0_FID),
                }
            )

    # null: M=0 ODE ⇒ isotropic
    sh0 = integrate_shells(r_f_grid, 0.0, D, eta_fid)
    fit0 = H_obs_from_radial_field(sh0["r_final"], sh0["v_final"], D, 0.1 * D, H_ref=H_f)
    null_ode = {
        "max_abs_frac_proxy": fit0["A_half_ptp"] + abs(fit0["A"]) + abs(fit0["Q_P2_moment"]),
        "A": fit0["A"],
        "Q": fit0["Q_P2_moment"],
        "half_ptp": fit0["A_half_ptp"],
        "pass": bool(fit0["A_half_ptp"] < 1e-3 and fit0["A"] < 1e-3),
    }

    # enhancement factor summary: integrated vs bolt-on α=1
    summary = {
        "eta_fid": eta_fid,
        "alpha_frozen_max": alpha_eff_frozen(0.0),
        "alpha_drag_max": alpha_eff_drag(0.0),
        "alpha_frozen_at_eta": a_fr,
        "alpha_drag_at_eta": a_dr,
        "kappa_ratio_frozen_over_bolt": 1.0 / a_fr,  # since A∝ακ, κ_needed ∝ 1/α
        "kappa_ratio_drag_over_bolt": 1.0 / a_dr,
        "conclusion_O1": (
            "Time integration on homologous ST changes α by O(1) only "
            f"(frozen α≈{a_fr:.2f}, drag α≈{a_dr:.2f} vs bolt-on α=1). "
            "It does NOT remove the need for κ~0.01–0.1 to get 5–15% sky variation."
        ),
    }

    return {
        "analytic_alpha_eff": analytic,
        "homology_null_ode": homology_null,
        "null_ode_M0": null_ode,
        "required_kappa_effective_alpha": required,
        "ode_scan": ode_rows,
        "ode_required": ode_req,
        "physical_masses": phys,
        "summary": summary,
        "H0": H0_FID,
        "D_H_Mpc": hubble_distance_mpc(H0_FID),
    }


if __name__ == "__main__":
    res = run_integrated_analysis()
    out = Path("/workspace/cosmology/integrated_results.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=2)
    print("Homology null:", res["homology_null_ode"])
    print("ODE M=0 null:", res["null_ode_M0"])
    print("Summary:", res["summary"])
    print("\nRequired κ (half-ptp metric), d/D=0.1:")
    for r in res["required_kappa_effective_alpha"]:
        if r["d_over_D"] == 0.1:
            print(
                f"  A~{r['target_half_ptp']}: bolt κ={r['kappa_bolt_alpha1']:.4f} "
                f"frozen κ={r['kappa_frozen_int']:.4f} drag κ={r['kappa_drag_int']:.4f}"
            )
    print("\nODE required:")
    for r in res["ode_required"]:
        print(f"  d/D={r['d_over_D']} T={r['target_half_ptp']}: κ_ode={r['kappa_ode']} ({r['note']})")
    print("\nPhysical (drag_int, d/D=0.1):")
    for p in res["physical_masses"]:
        if p["model"] == "drag_int":
            print(f"  target={p['target']} κ={p['kappa']:.4f}")
            for m in p["masses"]:
                print(
                    f"    D={m['D_Gpc']} Gpc: M={m['M_Msun']:.3e} Msun = "
                    f"{m['M_over_Mcrit']:.3f}×Mcrit(D) = {m['M_over_Mcrit_DH']:.3e}×Mcrit(DH)"
                )

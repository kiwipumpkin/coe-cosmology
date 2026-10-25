#!/usr/bin/env python3
"""
Dipole fraction vs d/D — observer-centric Migkas-style fits.

Addresses Chris pushback:
- Full-sky Legendre C1≈0 does NOT imply a Migkas A cosφ fit on mixed/partial
  sky returns A≈0; quadrupole-rich fields can leak into fitted dipoles.
- Measure C1,C2,C3 (orthogonal) AND unconstrained A cosφ (Migkas template)
  with χ² residual to show when the dipole description is good/bad.
- Hemisphere / cone selection: quantify fake dipole.

κ = κ_natural = 0.15, α = 1 (bolt-on); also α_eff_drag ≈ 2/3 for reference.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from forward_model import ModelParams, H_obs_map, make_latlon_grid
from integrated_gravity import alpha_eff_drag

OUT = Path("/workspace/cosmology")
KAPPA = 0.15
ALPHA = 1.0
DD_LIST = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
ETA_DRAG = 1e-4  # early start → α_eff → 2/3


def legendre_P(ell: int, mu: np.ndarray) -> np.ndarray:
    """Legendre P_ell(mu) via recurrence."""
    mu = np.asarray(mu, dtype=float)
    if ell == 0:
        return np.ones_like(mu)
    if ell == 1:
        return mu.copy()
    Pnm2 = np.ones_like(mu)
    Pnm1 = mu.copy()
    for n in range(2, ell + 1):
        Pn = ((2 * n - 1) * mu * Pnm1 - (n - 1) * Pnm2) / n
        Pnm2, Pnm1 = Pnm1, Pn
    return Pnm1


def frac_vs_mu(d_over_D: float, kappa: float, alpha: float, n_mu: int = 512):
    """Exact axisymmetric H_obs/H-1 on μ = n·D̂ grid (analytic geometry)."""
    D, H, G = 1.0, 1.0, 1.0
    d = d_over_D * D
    p = ModelParams(D=D, M=kappa, H=H, alpha=alpha, G=G)
    # μ from -1+eps to 1-eps; avoid CoE singularity when d≈D and μ=-1
    mu = np.linspace(-1.0 + 1e-6, 1.0 - 1e-6, n_mu)
    # n_hat in xz plane: n = (μ, 0, sqrt(1-μ²)) — or just (μ,0,0) wait need |n|=1
    # Use n = (μ, 0, sqrt(1-μ²))
    s = np.sqrt(np.maximum(1.0 - mu**2, 0.0))
    n_hat = np.stack([mu, np.zeros_like(mu), s], axis=-1)
    _, _, frac = H_obs_map(n_hat, d, p)
    return mu, frac


def legendre_coeffs(mu: np.ndarray, frac: np.ndarray, ell_max: int = 4) -> dict:
    """
    Project onto Legendre: frac(μ) = Σ C_ell P_ell(μ)
    with C_ell = (2ell+1)/2 ∫_{-1}^{1} frac(μ) P_ell(μ) dμ
    Use Gauss–Legendre quadrature on resampled interpolant.
    """
    # interpolate to Gauss nodes
    x, w = leggauss(256)
    # clip domain of mu
    m0, m1 = float(mu.min()), float(mu.max())
    # frac may not cover full [-1,1] exactly; interpolate within range, nan outside
    f_at = np.interp(x, mu, frac)
    # if endpoints not ±1, still integrate over available support ≈ [-1,1]
    coeffs = {}
    for ell in range(ell_max + 1):
        P = legendre_P(ell, x)
        C = (2 * ell + 1) / 2.0 * np.sum(w * f_at * P)
        coeffs[f"C{ell}"] = float(C)
    return coeffs


def fourier_equatorial(d_over_D: float, kappa: float, alpha: float, n_phi: int = 720):
    """H_obs(φ) on equator (lat=0): φ=0 toward CoE (+x)."""
    D, H, G = 1.0, 1.0, 1.0
    d = d_over_D * D
    p = ModelParams(D=D, M=kappa, H=H, alpha=alpha, G=G)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    n_hat = np.stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)], axis=-1)
    _, _, frac = H_obs_map(n_hat, d, p)
    # Discrete Fourier: a_k = (2/N) Σ f cos(kφ) for k≥1; a0 = mean
    # For real even function about φ=0 (axisym): only cos terms.
    N = n_phi
    a0 = float(np.mean(frac))
    cos_amps = {}
    sin_amps = {}
    for k in range(1, 5):
        ck = 2.0 / N * np.sum(frac * np.cos(k * phi))
        sk = 2.0 / N * np.sum(frac * np.sin(k * phi))
        cos_amps[f"a{k}"] = float(ck)
        sin_amps[f"b{k}"] = float(sk)
    return phi, frac, {"a0": a0, **cos_amps, **sin_amps}


def fit_pure_dipole_template(frac: np.ndarray, mu: np.ndarray, weights: np.ndarray | None = None):
    """
    Migkas-style: frac ≈ A * μ  (pure dipole, NO monopole in template),
    and also frac ≈ M + A * μ (dipole+monopole).
    Report A, χ²/ν, R², residual RMS, variance explained.
    """
    if weights is None:
        weights = np.ones_like(frac)
    w = weights / np.sum(weights)

    def wls(cols):
        X = np.column_stack(cols)
        sw = np.sqrt(w)
        coef, *_ = np.linalg.lstsq(X * sw[:, None], frac * sw, rcond=None)
        pred = X @ coef
        resid = frac - pred
        ss_res = float(np.sum(w * resid**2))
        ss_tot = float(np.sum(w * (frac - np.average(frac, weights=w)) ** 2))
        chi2 = float(np.sum(w * resid**2) * len(frac))  # scale-free proxy
        # Proper reduced χ² proxy: RMS_res / RMS_raw
        rms_res = float(np.sqrt(np.average(resid**2, weights=w)))
        rms_raw = float(np.sqrt(np.average(frac**2, weights=w)))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        return coef, pred, {
            "ss_res": ss_res,
            "ss_tot": ss_tot,
            "R2": r2,
            "rms_residual": rms_res,
            "rms_raw": rms_raw,
            "rms_ratio": rms_res / rms_raw if rms_raw > 0 else None,
            "chi2_proxy": ss_res / (ss_tot + 1e-30),  # 0=perfect, 1=no better than mean
        }

    # pure A μ (no monopole)
    coef1, pred1, m1 = wls([mu])
    # M + A μ
    coef2, pred2, m2 = wls([np.ones_like(mu), mu])
    # M + Q P2 (honest quadrupole)
    P2 = 0.5 * (3.0 * mu**2 - 1.0)
    coef3, pred3, m3 = wls([np.ones_like(mu), P2])
    # M + A μ + Q P2
    coef4, pred4, m4 = wls([np.ones_like(mu), mu, P2])

    return {
        "pure_dipole_A_mu": {
            "A": float(coef1[0]),
            **m1,
        },
        "monopole_plus_dipole": {
            "M": float(coef2[0]),
            "A": float(coef2[1]),
            **m2,
        },
        "monopole_plus_quad": {
            "M": float(coef3[0]),
            "Q": float(coef3[1]),
            **m3,
        },
        "monopole_dipole_quad": {
            "M": float(coef4[0]),
            "A": float(coef4[1]),
            "Q": float(coef4[2]),
            **m4,
        },
    }


def full_sky_analysis(d_over_D: float, kappa: float, alpha: float):
    mu, frac = frac_vs_mu(d_over_D, kappa, alpha, n_mu=1024)
    # Legendre with uniform dμ weight (axisymmetric ring average already 1D)
    C = legendre_coeffs(mu, frac, ell_max=4)
    fits = fit_pure_dipole_template(frac, mu, weights=np.ones_like(mu))
    phi, frac_eq, fourier = fourier_equatorial(d_over_D, kappa, alpha)
    # Migkas fit on equator: A cosφ
    fits_eq = fit_pure_dipole_template(frac_eq, np.cos(phi), weights=np.ones_like(phi))
    half_ptp = 0.5 * (float(frac.max()) - float(frac.min()))
    return {
        "d_over_D": d_over_D,
        "kappa": kappa,
        "alpha": alpha,
        "legendre": C,
        "C1_over_C2": abs(C["C1"] / C["C2"]) if abs(C["C2"]) > 1e-30 else None,
        "half_ptp": half_ptp,
        "fits_mu": fits,
        "fourier_eq": fourier,
        "fits_equator": fits_eq,
        "frac_mu_sample": {
            "mu": mu[::32].tolist(),
            "frac": frac[::32].tolist(),
        },
    }


def hemisphere_cone_analysis(d_over_D: float, kappa: float, alpha: float):
    """
    Full 3D sky grid; fit Migkas A (n·D̂) on:
      - full sky
      - CoE-facing hemisphere (μ>0)
      - anti-CoE hemisphere (μ<0)
      - cones: μ > cos(θ) for θ = 60°, 45°, 30° half-opening from CoE axis
    """
    D, H, G = 1.0, 1.0, 1.0
    d = d_over_D * D
    p = ModelParams(D=D, M=kappa, H=H, alpha=alpha, G=G)
    n_hat, lon, lat, w = make_latlon_grid(120, 60)
    _, _, frac = H_obs_map(n_hat, d, p)
    mu = n_hat[:, 0]  # CoE along +x

    selections = {
        "full_sky": np.ones_like(mu, dtype=bool),
        "CoE_facing_hemisphere_mu_gt_0": mu > 0,
        "anti_CoE_hemisphere_mu_lt_0": mu < 0,
        "cone_60deg_from_CoE": mu > np.cos(np.deg2rad(60)),
        "cone_45deg_from_CoE": mu > np.cos(np.deg2rad(45)),
        "cone_30deg_from_CoE": mu > np.cos(np.deg2rad(30)),
        "anti_cone_60deg": mu < -np.cos(np.deg2rad(60)),
    }

    out = {}
    for name, mask in selections.items():
        if mask.sum() < 10:
            out[name] = {"error": "too few pixels"}
            continue
        ww = w[mask].copy()
        ww = ww / ww.sum()
        fits = fit_pure_dipole_template(frac[mask], mu[mask], weights=ww)
        # also report dipole fraction: |A| / (|A|+|Q|) from joint fit if available
        mdq = fits["monopole_dipole_quad"]
        mpd = fits["monopole_plus_dipole"]
        mpq = fits["monopole_plus_quad"]
        A = mpd["A"]
        Q = mpq["Q"]
        out[name] = {
            "n_pix": int(mask.sum()),
            "sky_frac_weight": float(w[mask].sum()),
            "Migkas_A_from_M_plus_A_mu": A,
            "R2_dipole": mpd["R2"],
            "chi2_proxy_dipole": mpd["chi2_proxy"],
            "rms_ratio_dipole": mpd["rms_ratio"],
            "Q_from_M_plus_P2": Q,
            "R2_quad": mpq["R2"],
            "chi2_proxy_quad": mpq["chi2_proxy"],
            "joint_A": mdq["A"],
            "joint_Q": mdq["Q"],
            "R2_joint": mdq["R2"],
            "dipole_fraction_AQ": abs(A) / (abs(A) + abs(Q)) if (abs(A) + abs(Q)) > 0 else None,
            "fake_dipole_vs_Migkas_scale": abs(A) / 0.10,  # relative to A~0.1 target
        }
    return {"d_over_D": d_over_D, "kappa": kappa, "alpha": alpha, "selections": out}


def classify_migkas_fit(row: dict) -> str:
    """When is pure-dipole a good description?"""
    # Use equator M+A fit R2 and |C1|/|C2|
    r2 = row["fits_equator"]["monopole_plus_dipole"]["R2"]
    ratio = row["C1_over_C2"]
    if ratio is None:
        return "undefined"
    if ratio > 0.5 and r2 > 0.85:
        return "GOOD — dipole-dominated; Migkas A cosφ adequate"
    if ratio > 0.2 and r2 > 0.7:
        return "MARGINAL — mixed; dipole template incomplete"
    if r2 < 0.3 and ratio < 0.05:
        return "BAD — quadrupole-dominated; Migkas dipole mis-describes sky"
    return f"POOR — R2_dip={r2:.2f}, |C1/C2|={ratio:.3g}"


def main():
    results = {
        "kappa_natural": KAPPA,
        "alpha_bolt_on": ALPHA,
        "alpha_eff_drag_eta": ETA_DRAG,
        "alpha_eff_drag": float(alpha_eff_drag(ETA_DRAG)),
        "note": (
            "C_ell are full-sky Legendre coeffs of frac(μ)=H_obs/H-1. "
            "Orthogonality ⇒ C1≈0 at small d/D for point-mass tidal field. "
            "Migkas A cosφ is a SEPARATE fitted amplitude — can be nonzero on "
            "partial sky even when C1_full≈0. Do not equate 'quadrupole' with "
            "'no fitted dipole on any mask'."
        ),
        "full_sky_bolt_on": [],
        "full_sky_drag": [],
        "hemisphere_bolt_on": [],
    }

    a_drag = float(alpha_eff_drag(ETA_DRAG))

    for dD in DD_LIST:
        row = full_sky_analysis(dD, KAPPA, ALPHA)
        row["migkas_fit_quality"] = classify_migkas_fit(row)
        # drop bulky samples from summary tables but keep in JSON lightly
        results["full_sky_bolt_on"].append(row)

        row_d = full_sky_analysis(dD, KAPPA, a_drag)
        row_d["migkas_fit_quality"] = classify_migkas_fit(row_d)
        results["full_sky_drag"].append(row_d)

        hemi = hemisphere_cone_analysis(dD, KAPPA, ALPHA)
        results["hemisphere_bolt_on"].append(hemi)

    # Compact tables for report
    table_C = []
    for row in results["full_sky_bolt_on"]:
        C = row["legendre"]
        fe = row["fits_equator"]["monopole_plus_dipole"]
        fq = row["fits_equator"]["monopole_plus_quad"]
        table_C.append(
            {
                "d_over_D": row["d_over_D"],
                "C0": C["C0"],
                "C1": C["C1"],
                "C2": C["C2"],
                "C3": C["C3"],
                "C4": C["C4"],
                "C1_over_C2": row["C1_over_C2"],
                "Migkas_A_full_eq": fe["A"],
                "R2_dipole": fe["R2"],
                "chi2_proxy_dipole": fe["chi2_proxy"],
                "Q_full_eq": fq["Q"],
                "R2_quad": fq["R2"],
                "half_ptp": row["half_ptp"],
                "quality": row["migkas_fit_quality"],
                "fourier_a1": row["fourier_eq"]["a1"],
                "fourier_a2": row["fourier_eq"]["a2"],
            }
        )
    results["table_C1_C2"] = table_C

    table_hemi = []
    for hemi in results["hemisphere_bolt_on"]:
        dD = hemi["d_over_D"]
        for sel, v in hemi["selections"].items():
            if "error" in v:
                continue
            table_hemi.append(
                {
                    "d_over_D": dD,
                    "selection": sel,
                    "A_Migkas": v["Migkas_A_from_M_plus_A_mu"],
                    "R2_dipole": v["R2_dipole"],
                    "Q": v["Q_from_M_plus_P2"],
                    "R2_quad": v["R2_quad"],
                    "dipole_fraction_AQ": v["dipole_fraction_AQ"],
                    "A_over_0.1": v["fake_dipole_vs_Migkas_scale"],
                }
            )
    results["table_hemisphere"] = table_hemi

    # Strip heavy samples for JSON size
    for block in ("full_sky_bolt_on", "full_sky_drag"):
        for row in results[block]:
            row.pop("frac_mu_sample", None)

    out_path = OUT / "dipole_fraction_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    # Print summary
    print("=== Full-sky Legendre + Migkas A cosφ (κ=0.15, α=1) ===")
    print(
        f"{'d/D':>6} {'C1':>10} {'C2':>10} {'C3':>10} {'|C1/C2|':>10} "
        f"{'A_Mig':>10} {'R2_dip':>8} {'R2_Q':>8}  quality"
    )
    for t in table_C:
        print(
            f"{t['d_over_D']:6.2f} {t['C1']:10.3e} {t['C2']:10.4f} {t['C3']:10.3e} "
            f"{t['C1_over_C2']:10.3e} {t['Migkas_A_full_eq']:10.3e} "
            f"{t['R2_dipole']:8.3f} {t['R2_quad']:8.3f}  {t['quality'][:40]}"
        )

    print("\n=== Hemisphere / cone fake dipole (α=1) ===")
    print(f"{'d/D':>6} {'selection':<32} {'A':>10} {'R2_dip':>8} {'Q':>10} {'R2_Q':>8}")
    for t in table_hemi:
        if t["selection"] in (
            "full_sky",
            "CoE_facing_hemisphere_mu_gt_0",
            "anti_CoE_hemisphere_mu_lt_0",
            "cone_60deg_from_CoE",
        ):
            print(
                f"{t['d_over_D']:6.2f} {t['selection']:<32} {t['A_Migkas']:10.4f} "
                f"{t['R2_dipole']:8.3f} {t['Q']:10.4f} {t['R2_quad']:8.3f}"
            )

    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()

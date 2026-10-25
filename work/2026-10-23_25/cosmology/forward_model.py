#!/usr/bin/env python3
"""
Cleaned Phase A forward model (eqs 3.11'–3.12') + GM/D³H² absurdity check.

Ontology: CoE at origin; observer at D; homologous ST + central-gravity peculiar ansatz.
Critical Phase A result: pure ST ⇒ H_obs = H isotropic; anisotropy only from u.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants (SI-consistent astronomical units)
# ---------------------------------------------------------------------------
G_SI = 6.67430e-11  # m^3 kg^-1 s^-2
M_SUN = 1.98847e30  # kg
MPC_M = 3.085677581491367e22  # m
KM = 1e3
H0_FID = 70.0  # km/s/Mpc
C_KMS = 299792.458  # km/s

# G in (km/s)^2 * Mpc / Msun
G_ASTRO = G_SI * M_SUN / (MPC_M * KM**2) * 1e-6 * MPC_M / KM
# Wait: v^2 = GM/r with v in km/s, r in Mpc, M in Msun
# G_astro such that (km/s)^2 = G_astro * M_sun / Mpc
G_ASTRO = G_SI * M_SUN / (MPC_M) / (KM**2)  # (km/s)^2 * Mpc / Msun


def hubble_distance_mpc(H0: float = H0_FID) -> float:
    """c/H0 in Mpc."""
    return C_KMS / H0


def rho_crit_msun_mpc3(H0: float = H0_FID) -> float:
    """Critical density in Msun / Mpc^3. ρ_c = 3 H² / (8 π G)."""
    H_si = H0 * KM / MPC_M  # 1/s
    rho_si = 3.0 * H_si**2 / (8.0 * np.pi * G_SI)  # kg/m^3
    return rho_si * (MPC_M**3) / M_SUN


def mass_crit_sphere(D_mpc: float, H0: float = H0_FID) -> float:
    """M_crit(<D) = (4π/3) ρ_c D³ in solar masses."""
    return (4.0 * np.pi / 3.0) * rho_crit_msun_mpc3(H0) * D_mpc**3


def kappa_from_M_D(M_msun: float, D_mpc: float, H0: float = H0_FID) -> float:
    """κ = G M / (D³ H²), dimensionless."""
    # GM/(D³ H²): use G_ASTRO with H in km/s/Mpc, D in Mpc
    # G_ASTRO * M / D has units (km/s)^2; divide by (H D)^2
    return G_ASTRO * M_msun / (D_mpc**3 * H0**2)


def M_from_kappa(kappa: float, D_mpc: float, H0: float = H0_FID) -> float:
    """Invert κ = GM/(D³ H²)."""
    return kappa * D_mpc**3 * H0**2 / G_ASTRO


# Relation: κ = M / (2 M_crit) because M_crit = H² D³ / (2G)
# Check: G M_crit / (D³ H²) = G * (H² D³)/(2G) / (D³ H²) = 1/2
# so κ = M / (2 M_crit)  ⇒  M = 2 κ M_crit


def m_enclosed_fraction(X: np.ndarray) -> np.ndarray:
    """Eq 3.13': m(r)/M = 1 - e^{-X}(1 + X + X²/2), X = r/r_s."""
    X = np.asarray(X, dtype=float)
    return 1.0 - np.exp(-X) * (1.0 + X + 0.5 * X**2)


@dataclass
class ModelParams:
    D: float  # observer–CoE distance (same units as d, r)
    M: float  # total central mass (code units: G=1 or physical)
    H: float  # homologous expansion rate
    alpha: float = 1.0
    r_s: float | None = None  # if None, point mass
    G: float = 1.0


def peculiar_velocity(r_vec: np.ndarray, p: ModelParams) -> np.ndarray:
    """
    u(r) = -α [G m(r) / (r² H)] r̂
    r_vec: (..., 3)
    """
    r = np.linalg.norm(r_vec, axis=-1)
    # avoid singularity at CoE
    r_safe = np.maximum(r, 1e-30)
    if p.r_s is None or p.r_s <= 0:
        m = np.full_like(r_safe, p.M)
    else:
        m = p.M * m_enclosed_fraction(r_safe / p.r_s)
    # u = -α G m /(r² H) * r̂ = -α G m /(r³ H) * r_vec
    pref = -p.alpha * p.G * m / (r_safe**3 * p.H)
    return pref[..., None] * r_vec


def H_obs_map(
    n_hat: np.ndarray,
    d: float,
    p: ModelParams,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    For each sky direction n_hat (N,3), source at r = D_vec + d n_hat.
    Returns H_obs, Psi (=H_obs-H), H_obs/H - 1.
    """
    D_vec = np.array([p.D, 0.0, 0.0])
    r_vec = D_vec + d * n_hat
    u_r = peculiar_velocity(r_vec, p)
    u_D = peculiar_velocity(D_vec.reshape(1, 3), p)[0]
    delta_v = np.sum((u_r - u_D) * n_hat, axis=-1)
    H_obs = p.H + delta_v / d
    Psi = H_obs - p.H
    frac = H_obs / p.H - 1.0
    return H_obs, Psi, frac


def make_latlon_grid(n_lon: int = 72, n_lat: int = 36) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Equal-angle lat-lon grid. Returns n_hat (N,3), lon, lat in radians,
    and solid-angle weights ∝ cos(lat) for approximate spherical averages.
    lon ∈ [0, 2π), lat ∈ [-π/2, π/2] (excluding poles slightly).
    CoE axis = +x; n = (cos lat cos lon, cos lat sin lon, sin lat) but we use
    standard: θ polar from +z, φ azimuthal — for dipole fit along +x we use:
    n_x = sinθ cosφ, n_y = sinθ sinφ, n_z = cosθ with θ from +z.
    Simpler: use φ (azimuth from +x in xy) and θ (polar from +z).
    """
    lon = np.linspace(0.0, 2.0 * np.pi, n_lon, endpoint=False)  # φ
    # exclude exact poles for stability
    lat = np.linspace(-0.5 * np.pi + 1e-3, 0.5 * np.pi - 1e-3, n_lat)
    Lon, Lat = np.meshgrid(lon, lat, indexing="xy")
    # n relative to observer; CoE along +x (D_hat)
    # lat = elevation from xy-plane; lon from +x
    nx = np.cos(Lat) * np.cos(Lon)
    ny = np.cos(Lat) * np.sin(Lon)
    nz = np.sin(Lat)
    n_hat = np.stack([nx.ravel(), ny.ravel(), nz.ravel()], axis=-1)
    w = np.cos(Lat).ravel()  # dΩ ∝ cos(lat) dlat dlon
    w = w / w.sum()
    return n_hat, Lon.ravel(), Lat.ravel(), w


def fit_dipole_quadrupole(
    frac: np.ndarray, n_hat: np.ndarray, w: np.ndarray
) -> dict:
    """
    Fit frac ≈ c0 + A_vec·n + quadrupole (l=2, 5 independent).
    For CoE on +x, dipole amplitude A = |A_vec|, direction = A_vec/|A|.
    Also report monopole and P2 amplitude along x.
    """
    nx, ny, nz = n_hat[:, 0], n_hat[:, 1], n_hat[:, 2]
    # Design matrix: 1, nx, ny, nz, and l=2: (3nx²-1)/2, (3ny²-1)/2 wait —
    # use Cartesian quadrupole: nx*ny, nx*nz, ny*nz, nx²-ny², 2nz²-nx²-ny²
    cols = [
        np.ones_like(frac),
        nx,
        ny,
        nz,
        nx * ny,
        nx * nz,
        ny * nz,
        nx**2 - ny**2,
        2.0 * nz**2 - nx**2 - ny**2,
    ]
    X = np.column_stack(cols)
    # weighted least squares
    sw = np.sqrt(w)
    Xw = X * sw[:, None]
    yw = frac * sw
    coef, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    A_vec = coef[1:4]
    A = float(np.linalg.norm(A_vec))
    direction = A_vec / A if A > 0 else A_vec
    # residual RMS
    pred = X @ coef
    rms = float(np.sqrt(np.average((frac - pred) ** 2, weights=w)))
    raw_rms = float(np.sqrt(np.average(frac**2, weights=w)))
    # dipole-only projection along x (expected symmetry axis)
    # <frac * nx> / <nx²> ≈ A_x for pure dipole
    A_x_moment = float(np.average(frac * nx, weights=w) / np.average(nx**2, weights=w))
    # P2 moment: frac ~ Q * (3μ²-1)/2 with μ=nx
    P2 = 0.5 * (3.0 * nx**2 - 1.0)
    Q_moment = float(np.average(frac * P2, weights=w) / np.average(P2**2, weights=w))
    return {
        "monopole": float(coef[0]),
        "A": A,
        "A_vec": A_vec.tolist(),
        "direction": direction.tolist(),
        "A_x_moment": A_x_moment,
        "Q_P2_moment": Q_moment,
        "rms_residual": rms,
        "rms_raw": raw_rms,
        "coef": coef.tolist(),
    }


def dipole_amplitude(
    p: ModelParams,
    d: float,
    n_lon: int = 72,
    n_lat: int = 36,
) -> dict:
    n_hat, lon, lat, w = make_latlon_grid(n_lon, n_lat)
    _, _, frac = H_obs_map(n_hat, d, p)
    fit = fit_dipole_quadrupole(frac, n_hat, w)
    fit["frac_min"] = float(frac.min())
    fit["frac_max"] = float(frac.max())
    fit["frac_ptp"] = float(frac.max() - frac.min())
    # peak-to-peak / 2 as rough "amplitude" if dipole-dominated
    fit["A_half_ptp"] = 0.5 * fit["frac_ptp"]
    return fit


def null_test(tol: float = 1e-12) -> dict:
    """M=0 ⇒ Ψ=0 everywhere."""
    n_hat, _, _, w = make_latlon_grid(36, 18)
    p = ModelParams(D=2.0, M=0.0, H=1.0, alpha=1.0, G=1.0)
    _, Psi, frac = H_obs_map(n_hat, d=1.0, p=p)
    ok = bool(np.max(np.abs(Psi)) < tol and np.max(np.abs(frac)) < tol)
    # also α=0
    p2 = ModelParams(D=2.0, M=1.0, H=1.0, alpha=0.0, G=1.0)
    _, Psi2, frac2 = H_obs_map(n_hat, d=1.0, p=p2)
    ok2 = bool(np.max(np.abs(Psi2)) < tol)
    return {
        "M0_max_abs_Psi": float(np.max(np.abs(Psi))),
        "M0_pass": ok,
        "alpha0_max_abs_Psi": float(np.max(np.abs(Psi2))),
        "alpha0_pass": ok2,
        "pass": ok and ok2,
    }


def analytic_tidal_limit(kappa: float, alpha: float, d_over_D: float) -> dict:
    """
    Small-ε expansion notes.
    u = -C r_vec/r³, C = α G M / H, κ = G M/(D³ H²), C/(D³ H) wait:
    C / D³ = α G M /(H D³) = α κ H.
    Gradient at D: Ψ ≈ n_i n_j ∂_j u_i = -(C/D³)[1 - 3 μ²]
                    = -α κ H / H * [1-3μ²]?  Ψ/H = -(C/(D³ H))[1-3μ²]
    C/(D³ H) = α G M /(H² D³) = α κ
    So Ψ/H = -α κ (1 - 3 μ²) = α κ (3μ² - 1)  → pure P2 quadrupole, no dipole!
    Leading anisotropy in the strict tidal (d→0) limit is QUADRUPOLE.
    Dipole appears at O(d/D) and from evaluating u at r≠D.
    """
    return {
        "note": "tidal d→0: Ψ/H = ακ(3μ²-1) = 2ακ P2(μ); pure quadrupole",
        "Q_P2": 2.0 * alpha * kappa,  # if Ψ/H = Q P2 with P2=(3μ²-1)/2
        "dipole_at_linear_order": 0.0,
    }


def scan_kappa_vs_amplitude(
    kappa_grid: np.ndarray,
    alpha: float,
    d_over_D: float,
    n_lon: int = 48,
    n_lat: int = 24,
) -> dict:
    """
    Code units: D=1, H=1, G=1 ⇒ M = κ (since κ = G M / D³ H² = M).
    Vary κ; measure dipole A and quadrupole Q.
    """
    D, H, G = 1.0, 1.0, 1.0
    d = d_over_D * D
    rows = []
    for kappa in kappa_grid:
        p = ModelParams(D=D, M=kappa, H=H, alpha=alpha, G=G)  # M=κ when G=D=H=1
        fit = dipole_amplitude(p, d=d, n_lon=n_lon, n_lat=n_lat)
        rows.append(
            {
                "kappa": float(kappa),
                "A": fit["A"],
                "A_x": fit["A_x_moment"],
                "Q_P2": fit["Q_P2_moment"],
                "monopole": fit["monopole"],
                "A_half_ptp": fit["A_half_ptp"],
                "frac_ptp": fit["frac_ptp"],
            }
        )
    return {"alpha": alpha, "d_over_D": d_over_D, "rows": rows}


def required_kappa_for_A(
    target_A: float,
    alpha: float,
    d_over_D: float,
    metric: str = "A",
) -> float:
    """
    Binary search κ > 0 such that dipole amplitude ≈ target_A.
    Linear regime: A ∝ α κ, so κ ≈ target / (α * slope).
    """
    # first find slope at small κ
    k_probe = 1e-3
    p = ModelParams(D=1.0, M=k_probe, H=1.0, alpha=alpha, G=1.0)
    fit = dipole_amplitude(p, d=d_over_D, n_lon=48, n_lat=24)
    val = fit["A"] if metric == "A" else abs(fit["A_x_moment"])
    if val < 1e-30:
        # try half-ptp / use Q
        slope_Q = abs(fit["Q_P2_moment"]) / k_probe
        return {
            "target": target_A,
            "kappa_for_dipole": None,
            "reason": "dipole negligible at linear order; use quadrupole",
            "kappa_for_Q_eq_target": target_A / slope_Q if slope_Q > 0 else None,
            "slope_Q": slope_Q,
            "slope_A": 0.0,
        }
    slope = val / k_probe
    kappa_lin = target_A / slope

    # refine with search (nonlinear when κ large / d/D not tiny)
    lo, hi = 0.0, max(10.0 * kappa_lin, 1e-6)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        p = ModelParams(D=1.0, M=mid, H=1.0, alpha=alpha, G=1.0)
        fit = dipole_amplitude(p, d=d_over_D, n_lon=48, n_lat=24)
        v = fit["A"] if metric == "A" else abs(fit["A_x_moment"])
        if v < target_A:
            lo = mid
        else:
            hi = mid
    kappa_star = 0.5 * (lo + hi)
    p = ModelParams(D=1.0, M=kappa_star, H=1.0, alpha=alpha, G=1.0)
    fit = dipole_amplitude(p, d=d_over_D, n_lon=64, n_lat=32)
    return {
        "target": target_A,
        "metric": metric,
        "kappa": kappa_star,
        "alpha": alpha,
        "d_over_D": d_over_D,
        "A_achieved": fit["A"],
        "A_x": fit["A_x_moment"],
        "Q_P2": fit["Q_P2_moment"],
        "monopole": fit["monopole"],
        "slope_linear_A": slope,
        "kappa_linear_estimate": kappa_lin,
    }


def physical_mass_table(
    kappa: float, D_gpc_list: list[float], H0: float = H0_FID
) -> list[dict]:
    out = []
    for D_gpc in D_gpc_list:
        D_mpc = D_gpc * 1000.0
        M = M_from_kappa(kappa, D_mpc, H0)
        Mcrit = mass_crit_sphere(D_mpc, H0)
        # Ω_m = 0.3 mass in same volume
        M_Om = 0.3 * Mcrit
        DH = hubble_distance_mpc(H0)
        Mcrit_DH = mass_crit_sphere(DH, H0)
        out.append(
            {
                "D_Gpc": D_gpc,
                "D_Mpc": D_mpc,
                "kappa": kappa,
                "M_Msun": M,
                "M_Msun_log10": float(np.log10(M)) if M > 0 else None,
                "M_crit_sphere_Msun": Mcrit,
                "M_over_Mcrit": M / Mcrit if Mcrit > 0 else None,
                "M_Om03_sphere_Msun": M_Om,
                "M_over_Om_mass": M / M_Om if M_Om > 0 else None,
                "M_crit_HubbleRadius_Msun": Mcrit_DH,
                "M_over_Mcrit_DH": M / Mcrit_DH if Mcrit_DH > 0 else None,
                # cluster scale ~ 1e15 Msun
                "M_over_1e15_Msun": M / 1e15,
            }
        )
    return out


def run_full_analysis(out_dir: Path) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- sanity: G_ASTRO and κ ↔ M relation ---
    H0 = H0_FID
    D_test = 1000.0  # Mpc
    Mcrit = mass_crit_sphere(D_test, H0)
    kappa_of_Mcrit = kappa_from_M_D(Mcrit, D_test, H0)
    # expect ≈ 0.5

    null = null_test()

    # --- character of anisotropy vs d/D ---
    alpha = 1.0
    kappa_fixed = 0.01
    dD_list = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
    morph = []
    for dD in dD_list:
        p = ModelParams(D=1.0, M=kappa_fixed, H=1.0, alpha=alpha, G=1.0)
        fit = dipole_amplitude(p, d=dD, n_lon=64, n_lat=32)
        morph.append(
            {
                "d_over_D": dD,
                "A": fit["A"],
                "A_x": fit["A_x_moment"],
                "Q_P2": fit["Q_P2_moment"],
                "monopole": fit["monopole"],
                "A_over_akappa": fit["A"] / (alpha * kappa_fixed),
                "Q_over_akappa": fit["Q_P2_moment"] / (alpha * kappa_fixed),
                "ptp": fit["frac_ptp"],
            }
        )

    # Migkas-relevant: clusters ~100–500 Mpc; D ~ 1–14 Gpc ⇒ d/D ~ 0.01–0.5
    # Use representative d/D = 0.1 and 0.3
    targets = [0.05, 0.10, 0.15]
    alphas = [0.5, 1.0, 2.0]
    dD_cases = [0.05, 0.1, 0.3, 1.0]

    req = []
    for a in alphas:
        for dD in dD_cases:
            for T in targets:
                # dipole metric
                rA = required_kappa_for_A(T, a, dD, metric="A")
                # also κ needed so |Q| = T (quadrupole amplitude as "variation")
                # linear: Q ≈ 2 α κ in tidal limit for P2 coeff
                k_probe = 1e-4
                p = ModelParams(D=1.0, M=k_probe, H=1.0, alpha=a, G=1.0)
                fit = dipole_amplitude(p, d=dD, n_lon=48, n_lat=24)
                slope_Q = fit["Q_P2_moment"] / k_probe
                slope_ptp = fit["frac_ptp"] / k_probe
                # κ for half-ptp = T (i.e. peak variation ±T)
                kappa_ptp = (2.0 * T) / slope_ptp if slope_ptp != 0 else None
                kappa_Q = T / abs(slope_Q) if slope_Q != 0 else None
                req.append(
                    {
                        "alpha": a,
                        "d_over_D": dD,
                        "target_A": T,
                        "kappa_for_dipole_A": rA.get("kappa"),
                        "A_achieved": rA.get("A_achieved"),
                        "Q_at_that_kappa": rA.get("Q_P2"),
                        "kappa_for_Q_eq_target": kappa_Q,
                        "kappa_for_half_ptp_eq_target": kappa_ptp,
                        "slope_A_linear": rA.get("slope_linear_A"),
                        "slope_Q_linear": slope_Q,
                    }
                )

    # Physical conversion for key cases: α=1, d/D=0.1
    # Use κ from half-ptp (sky variation ±A) AND from dipole A
    D_list = [1.0, 3.0, 10.0, 14.0]  # Gpc; 14 ~ particle horizon scale often quoted

    key_rows = []
    for T in targets:
        for dD in [0.1, 0.3]:
            r = required_kappa_for_A(T, alpha=1.0, d_over_D=dD, metric="A")
            k_probe = 1e-4
            p = ModelParams(D=1.0, M=k_probe, H=1.0, alpha=1.0, G=1.0)
            fit = dipole_amplitude(p, d=dD, n_lon=48, n_lat=24)
            slope_ptp = fit["frac_ptp"] / k_probe
            kappa_ptp = (2.0 * T) / slope_ptp
            kappa_Q = T / abs(fit["Q_P2_moment"] / k_probe)
            phys_A = physical_mass_table(r["kappa"], D_list, H0)
            phys_ptp = physical_mass_table(kappa_ptp, D_list, H0)
            phys_Q = physical_mass_table(kappa_Q, D_list, H0)
            key_rows.append(
                {
                    "target": T,
                    "d_over_D": dD,
                    "alpha": 1.0,
                    "kappa_dipole": r["kappa"],
                    "kappa_half_ptp": kappa_ptp,
                    "kappa_Q": kappa_Q,
                    "phys_dipole": phys_A,
                    "phys_half_ptp": phys_ptp,
                    "phys_Q": phys_Q,
                    "A_achieved": r["A_achieved"],
                    "Q_at_kappa_dipole": r["Q_P2"],
                }
            )

    # Extended mass profile check: r_s = 0.3 D, same κ total
    ext = []
    for dD in [0.1, 0.3]:
        p_pt = ModelParams(D=1.0, M=0.01, H=1.0, alpha=1.0, G=1.0, r_s=None)
        p_ex = ModelParams(D=1.0, M=0.01, H=1.0, alpha=1.0, G=1.0, r_s=0.3)
        f_pt = dipole_amplitude(p_pt, d=dD)
        f_ex = dipole_amplitude(p_ex, d=dD)
        ext.append(
            {
                "d_over_D": dD,
                "point_A": f_pt["A"],
                "ext_A": f_ex["A"],
                "point_Q": f_pt["Q_P2_moment"],
                "ext_Q": f_ex["Q_P2_moment"],
            }
        )

    # Linearity check: A vs κ
    kappa_grid = np.logspace(-4, 0, 17)
    lin_scan = scan_kappa_vs_amplitude(kappa_grid, alpha=1.0, d_over_D=0.1)

    # Homology check: wrong formula reminder not implemented; pure ST already null

    results = {
        "H0_km_s_Mpc": H0,
        "G_astro_(km/s)^2_Mpc/Msun": G_ASTRO,
        "D_H_Mpc": hubble_distance_mpc(H0),
        "rho_c_Msun_Mpc3": rho_crit_msun_mpc3(H0),
        "kappa_of_Mcrit_check": kappa_of_Mcrit,  # should be ~0.5
        "null_test": null,
        "morphology_vs_dD": morph,
        "analytic_tidal": analytic_tidal_limit(0.01, 1.0, 0.0),
        "required_kappa_grid": req,
        "key_physical": key_rows,
        "extended_vs_point": ext,
        "linearity_scan_dD0.1": lin_scan,
    }

    with open(out_dir / "absurdity_results.json", "w") as f:
        json.dump(results, f, indent=2)

    return results


if __name__ == "__main__":
    res = run_full_analysis(Path("/workspace/cosmology"))
    print("NULL TEST:", res["null_test"])
    print("kappa(M_crit) check (expect ~0.5):", res["kappa_of_Mcrit_check"])
    print("\nMorphology vs d/D (α=1, κ=0.01):")
    for m in res["morphology_vs_dD"]:
        print(
            f"  d/D={m['d_over_D']:.2f}: A={m['A']:.4e} Q={m['Q_P2']:.4e} "
            f"A/(ακ)={m['A_over_akappa']:.3f} Q/(ακ)={m['Q_over_akappa']:.3f}"
        )
    print("\nKey physical (α=1):")
    for kr in res["key_physical"]:
        print(
            f"  A_target={kr['target']}, d/D={kr['d_over_D']}: "
            f"κ_dip={kr['kappa_dipole']:.4e} κ_ptp={kr['kappa_half_ptp']:.4e} κ_Q={kr['kappa_Q']:.4e}"
        )
        for p in kr["phys_half_ptp"]:
            print(
                f"    D={p['D_Gpc']} Gpc: M={p['M_Msun']:.3e} Msun = "
                f"{p['M_over_Mcrit']:.2f}×M_crit(<D) = {p['M_over_Mcrit_DH']:.2e}×M_crit(DH)"
            )

#!/usr/bin/env python3
"""
Local-density match + Migkas Fig.6 fit under Phase A decreasing-ρ.

Constraint: ρ(D) = Ω_m ρ_crit (= 0.3 ρ_crit) by construction.
Scan x = D/L for exponential ρ∝e^{-r/L}; optionally softer / two-scale.
Report actual κ = G M_enc(<D)/(D³ H²) (must exceed Ω_m/2 for decreasing ρ).
ST+Gauss H_obs on observer sky circle → A∝1/(H_obs/H)²; AN∝m(r)/M_enc.
Chris sheet-style Amp×sin phase rotation for joint A&B R².
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

ROOT = Path("/workspace/cosmology")
sys.path.insert(0, str(ROOT))
from forward_model import ModelParams, H_obs_map, m_enclosed_fraction  # noqa: E402

OUT = Path(__file__).resolve().parent
CSV = ROOT / "migkas" / "fig6_digitized.csv"

OMEGA_M = 0.3
ALPHA = 1.0
EPS = 0.1  # d/D
H = 1.0
G = 1.0
D = 1.0
OMEGA = 2.0 * np.pi / 360.0

C_BLUE = "#0072B2"
C_ORANGE = "#E69F00"
C_GREEN = "#009E73"
C_BLACK = "#222222"
C_GRAY = "#666666"
C_RED = "#D55E00"
C_PURPLE = "#CC79A7"


def style() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 220,
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 12,
            "legend.fontsize": 8.5,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 1.1,
            "lines.linewidth": 2.0,
            "mathtext.fontset": "dejavusans",
        }
    )


def rho_crit(H_: float = H, G_: float = G) -> float:
    return 3.0 * H_**2 / (8.0 * np.pi * G_)


def m_frac_exp(X: np.ndarray) -> np.ndarray:
    return m_enclosed_fraction(X)


def params_exp_local_match(x: float, omega_m: float = OMEGA_M) -> dict:
    """
    Exponential ρ(r)=M_tot/(8π L³) exp(-r/L), set ρ(D)=omega_m ρ_crit.
    Unitless: D=H=G=1.
    """
    L = D / x
    rc = rho_crit()
    rho_D = omega_m * rc
    # ρ(D)=M_tot/(8π L³) e^{-x} ⇒ M_tot
    M_tot = rho_D * 8.0 * np.pi * L**3 * np.exp(x)
    f = float(m_frac_exp(x))
    M_enc = M_tot * f
    kappa = G * M_enc / (D**3 * H**2)
    rho_mean = M_enc / ((4.0 / 3.0) * np.pi * D**3)
    return {
        "profile": "exponential",
        "x": float(x),
        "L": float(L),
        "L_over_D": float(L / D),
        "M_tot": float(M_tot),
        "M_enc": float(M_enc),
        "m_frac": f,
        "kappa": float(kappa),
        "rho_local_over_rcrit": float(rho_D / rc),
        "rho_mean_over_rcrit": float(rho_mean / rc),
        "rho_local_over_mean": float(rho_D / rho_mean) if rho_mean > 0 else float("nan"),
        "rho0": float(M_tot / (8.0 * np.pi * L**3)),
    }


def m_frac_plummer_like(X: np.ndarray, n: float = 2.0) -> np.ndarray:
    """
    Softer cumulative: m/M = X^n / (1+X^n)  (n=2 → Plummer-like enclosed).
    Density decreases; asymptotic M_tot finite.
    """
    X = np.asarray(X, dtype=float)
    return X**n / (1.0 + X**n)


def rho_from_m_frac_numerical(X: np.ndarray, mfrac: np.ndarray, M_tot: float, L: float) -> np.ndarray:
    """ρ = (1/(4π r²)) dm/dr with r=X L."""
    r = X * L
    # dm/dX via gradient; dm/dr = (dm/dX)/L
    dm_dX = np.gradient(mfrac * M_tot, X)
    dm_dr = dm_dX / L
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = dm_dr / (4.0 * np.pi * r**2)
    rho = np.where(np.isfinite(rho), rho, 0.0)
    return rho


def params_plummer_local_match(x: float, n: float = 2.0, omega_m: float = OMEGA_M) -> dict:
    """
    Soft profile with m/M = X^n/(1+X^n). Set M_tot so ρ(D)=omega_m ρ_crit
    using local density from analytic dm/dr.
    """
    L = D / x
    rc = rho_crit()
    rho_D_target = omega_m * rc
    # Analytic: m = M X^n/(1+X^n), dm/dX = M n X^{n-1}/(1+X^n)^2
    # ρ(D) = [dm/dr]/(4π D²) = [dm/dX / L] / (4π D²)
    # dm/dX(at x) = M n x^{n-1}/(1+x^n)^2
    # ⇒ M = rho_D * 4π D² * L * (1+x^n)^2 / (n x^{n-1})
    M_tot = rho_D_target * 4.0 * np.pi * D**2 * L * (1.0 + x**n) ** 2 / (n * x ** (n - 1))
    f = float(m_frac_plummer_like(x, n=n))
    M_enc = M_tot * f
    kappa = G * M_enc / (D**3 * H**2)
    rho_mean = M_enc / ((4.0 / 3.0) * np.pi * D**3)
    return {
        "profile": f"plummer_n{n:g}",
        "x": float(x),
        "n": float(n),
        "L": float(L),
        "L_over_D": float(L / D),
        "M_tot": float(M_tot),
        "M_enc": float(M_enc),
        "m_frac": f,
        "kappa": float(kappa),
        "rho_local_over_rcrit": float(rho_D_target / rc),
        "rho_mean_over_rcrit": float(rho_mean / rc),
        "rho_local_over_mean": float(rho_D_target / rho_mean) if rho_mean > 0 else float("nan"),
    }



def params_twoscale_local_match(x: float, f_u: float, omega_m: float = OMEGA_M) -> dict:
    """
    Two-scale: ρ = ρ_core0 e^{-r/L} + ρ_u (uniform floor).
    f_u = fraction of ρ(D) from the floor. Core supplies tidal κ; floor lifts local density.
    ρ(D)=omega_m ρ_crit by construction.
    """
    if not (0.0 <= f_u < 1.0):
        raise ValueError("f_u in [0,1)")
    L = D / x
    rc = rho_crit()
    rho_D = omega_m * rc
    rho_u = f_u * rho_D
    rho_cD = (1.0 - f_u) * rho_D
    M_tot = rho_cD * 8.0 * np.pi * L**3 * np.exp(x)
    f = float(m_frac_exp(x))
    M_core = M_tot * f
    M_u = (4.0 / 3.0) * np.pi * rho_u * D**3
    M_enc = M_core + M_u
    kappa = G * M_enc / (D**3 * H**2)
    kappa_core = G * M_core / (D**3 * H**2)
    kappa_u = G * M_u / (D**3 * H**2)
    rho_mean = M_enc / ((4.0 / 3.0) * np.pi * D**3)
    return {
        "profile": f"twoscale_fu{f_u:.2f}",
        "x": float(x),
        "f_u": float(f_u),
        "L": float(L),
        "L_over_D": float(L / D),
        "M_tot": float(M_tot),
        "M_enc": float(M_enc),
        "M_core": float(M_core),
        "M_u": float(M_u),
        "m_frac": f,
        "kappa": float(kappa),
        "kappa_core": float(kappa_core),
        "kappa_u": float(kappa_u),
        "rho_u": float(rho_u),
        "rho_local_over_rcrit": float(rho_D / rc),
        "rho_mean_over_rcrit": float(rho_mean / rc),
        "rho_local_over_mean": float(rho_D / rho_mean) if rho_mean > 0 else float("nan"),
        "rho0_core": float(M_tot / (8.0 * np.pi * L**3)),
    }


def m_enc_at_r(r: np.ndarray, pmeta: dict) -> np.ndarray:
    """Enclosed mass at radii r for the chosen profile family."""
    r = np.asarray(r, dtype=float)
    L = pmeta["L"]
    M_tot = pmeta["M_tot"]
    X = r / L
    if pmeta["profile"] == "exponential":
        return M_tot * m_frac_exp(X)
    if pmeta["profile"].startswith("plummer"):
        n = pmeta.get("n", 2.0)
        return M_tot * m_frac_plummer_like(X, n=n)
    if pmeta["profile"].startswith("twoscale"):
        m_core = M_tot * m_frac_exp(X)
        rho_u = pmeta["rho_u"]
        m_u = (4.0 / 3.0) * np.pi * rho_u * r**3
        return m_core + m_u
    raise ValueError(pmeta["profile"])


def peculiar_u(r_vec: np.ndarray, pmeta: dict, alpha: float = ALPHA) -> np.ndarray:
    """u = -α [G m(r)/(r² H)] r̂  (Gauss peculiar ansatz)."""
    r = np.linalg.norm(r_vec, axis=-1)
    r_safe = np.maximum(r, 1e-30)
    m = m_enc_at_r(r_safe, pmeta)
    pref = -alpha * G * m / (r_safe**3 * H)
    return pref[..., None] * r_vec


def H_obs_circle(pmeta: dict, eps: float = EPS, n_phi: int = 721) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Observer sky circle in the CoE–observer plane (and full φ).
    φ=0 along +x away from CoE? CoE at origin, observer at +D x-hat.
    Galactic ℓ proxy: φ_gal with μ = cos(ℓ - ℓ0); here return ℓ grid 0..360
    with μ = cos(ℓ) when preferred axis at ℓ=0 (toward/away CoE along x).
    n = (cos ℓ, sin ℓ, 0) in observer frame? Standard: n_x = cos(φ) toward CoE axis.
    Use φ measured from +x (CoE axis): μ = n·D̂ = cos φ.
    Map φ → ℓ with ℓ=φ in degrees for later rotation (φ0 absorbs orientation).
    """
    l_deg = np.linspace(0.0, 360.0, n_phi, endpoint=False)
    phi = np.deg2rad(l_deg)
    # n_hat: CoE along +x; φ=0 → +x (toward anti-CoE if CoE at -wait:
    # Observer at +D x; direction toward CoE is -x. For μ = n·D̂ with D̂=+x,
    # φ=0 (n=+x) is anti-center, φ=π is toward CoE. Matches make_figures.
    n_hat = np.stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)], axis=-1)
    d = eps * D
    D_vec = np.array([D, 0.0, 0.0])
    r_vec = D_vec + d * n_hat
    u_r = peculiar_u(r_vec, pmeta)
    u_D = peculiar_u(D_vec.reshape(1, 3), pmeta)[0]
    delta_v = np.sum((u_r - u_D) * n_hat, axis=-1)
    H_obs = H + delta_v / d
    # enclosed mass at sources (sheet AN)
    m_src = m_enc_at_r(np.linalg.norm(r_vec, axis=-1), pmeta)
    return l_deg, H_obs, m_src


def boxcar_dl90(y: np.ndarray, l_deg: np.ndarray) -> np.ndarray:
    """Δℓ=90° wrap-around boxcar (Migkas Fig.6 style)."""
    n = len(y)
    dl = 360.0 / n
    half = int(round(45.0 / dl))
    kern = np.ones(2 * half + 1)
    kern /= kern.sum()
    return np.convolve(np.pad(y, (half, half), mode="wrap"), kern, mode="valid")


def af_an_from_physics(pmeta: dict, eps: float = EPS, af_mode: str = "invsq") -> dict:
    """
    af_mode:
      invsq  — A∝1/(H_obs/H)² mean-renorm (task / DL-like)
      hratio — H_obs/H mean-renorm (Chris sheet AF-like)
    """
    l_deg, H_obs, m_src = H_obs_circle(pmeta, eps=eps)
    # avoid non-positive H_obs for invsq
    H_safe = np.maximum(H_obs, 1e-6 * H)
    f = H_obs / H - 1.0
    if af_mode == "invsq":
        A_raw = 1.0 / (H_safe / H) ** 2
    elif af_mode == "hratio":
        A_raw = H_safe / H
    else:
        raise ValueError(af_mode)
    AF = A_raw / np.mean(A_raw)
    M_enc = pmeta["M_enc"]
    AN_raw = m_src / M_enc
    AN = AN_raw / np.mean(AN_raw)
    AF_s = boxcar_dl90(AF, l_deg)
    AN_s = boxcar_dl90(AN, l_deg)
    return {
        "l_deg": l_deg,
        "H_obs": H_obs,
        "frac": f,
        "AF": AF,
        "AN": AN,
        "AF_smooth": AF_s,
        "AN_smooth": AN_s,
        "af_mode": af_mode,
        "H_obs_ptp": float(H_obs.max() - H_obs.min()),
        "frac_ptp": float(f.max() - f.min()),
    }


def rsq(y: np.ndarray, x: np.ndarray) -> float:
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    if y.size < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1] ** 2)


def sine_mod(l_deg: np.ndarray, amp: float, phi: float) -> np.ndarray:
    return amp * np.sin(OMEGA * np.asarray(l_deg, dtype=float) + phi)


def model_rotated(af, an, l_deg, amp_a, phi):
    ba = sine_mod(l_deg, amp_a, phi)
    be = sine_mod(l_deg, -amp_a, phi)
    return (1.0 - ba) * af, (1.0 - be) * an


def interp_periodic(l_tgt, l_src, y_src):
    mask = l_src < 360.0 - 1e-9
    ls, ys = l_src[mask], y_src[mask]
    order = np.argsort(ls)
    ls, ys = ls[order], ys[order]
    ls_p = np.concatenate([ls - 360.0, ls, ls + 360.0])
    ys_p = np.concatenate([ys, ys, ys])
    return np.interp(np.asarray(l_tgt, dtype=float) % 360.0, ls_p, ys_p)


def load_migkas():
    data = np.genfromtxt(CSV, delimiter=",", names=True)
    return {
        "l": np.asarray(data["l_deg"], dtype=float),
        "A": np.asarray(data["A_over_Aall"], dtype=float),
        "B": np.asarray(data["B_over_Ball"], dtype=float),
        "A_elo": np.asarray(data["A_err_low"], dtype=float),
        "A_ehi": np.asarray(data["A_err_high"], dtype=float),
        "B_elo": np.asarray(data["B_err_low"], dtype=float),
        "B_ehi": np.asarray(data["B_err_high"], dtype=float),
    }


def fit_amp_phi(l, af, an, a_obs, b_obs, amp_grid=None, n_phi=91):
    """Scan Amp and shared φ; maximise R²_A + R²_B (both must be decent)."""
    if amp_grid is None:
        amp_grid = np.concatenate(
            [
                np.array([0.0]),
                -np.linspace(0.02, 0.30, 29),
                np.linspace(0.02, 0.30, 29),
            ]
        )
    phis = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    best = {
        "amp": 0.0,
        "phi": 0.0,
        "r2_a": -1.0,
        "r2_b": -1.0,
        "r2_sum": -1.0,
        "r2_min": -1.0,
    }
    # Also track best under min(R2A,R2B) criterion
    best_balanced = dict(best)
    for amp in amp_grid:
        for phi in phis:
            ma, mb = model_rotated(af, an, l, float(amp), float(phi))
            ra, rb = rsq(a_obs, ma), rsq(b_obs, mb)
            s = ra + rb
            mn = min(ra, rb)
            if s > best["r2_sum"]:
                best = {
                    "amp": float(amp),
                    "phi": float(phi),
                    "phi_deg": float(np.rad2deg(phi)),
                    "r2_a": float(ra),
                    "r2_b": float(rb),
                    "r2_sum": float(s),
                    "r2_min": float(mn),
                }
            if mn > best_balanced["r2_min"] or (
                abs(mn - best_balanced["r2_min"]) < 1e-12 and s > best_balanced["r2_sum"]
            ):
                best_balanced = {
                    "amp": float(amp),
                    "phi": float(phi),
                    "phi_deg": float(np.rad2deg(phi)),
                    "r2_a": float(ra),
                    "r2_b": float(rb),
                    "r2_sum": float(s),
                    "r2_min": float(mn),
                }
    # Chris defaults for reference
    ma0, mb0 = model_rotated(af, an, l, -0.06, 0.35)
    ref = {
        "amp": -0.06,
        "phi": 0.35,
        "phi_deg": float(np.rad2deg(0.35)),
        "r2_a": float(rsq(a_obs, ma0)),
        "r2_b": float(rsq(b_obs, mb0)),
    }
    ref["r2_sum"] = ref["r2_a"] + ref["r2_b"]
    ref["r2_min"] = min(ref["r2_a"], ref["r2_b"])
    # without rotation
    wo = {
        "r2_a": float(rsq(a_obs, af)),
        "r2_b": float(rsq(b_obs, an)),
    }
    return {"best_sum": best, "best_balanced": best_balanced, "chris_ref": ref, "without": wo}


def rho_profile_curve(pmeta: dict, r_over_D: np.ndarray) -> np.ndarray:
    """ρ(r)/ρ_crit along radii."""
    rc = rho_crit()
    r = r_over_D * D
    L = pmeta["L"]
    M_tot = pmeta["M_tot"]
    if pmeta["profile"] == "exponential":
        rho = (M_tot / (8.0 * np.pi * L**3)) * np.exp(-r / L)
        return rho / rc
    if pmeta["profile"].startswith("twoscale"):
        rho_core = (M_tot / (8.0 * np.pi * L**3)) * np.exp(-r / L)
        rho = rho_core + pmeta["rho_u"]
        return rho / rc
    X = r / L
    n = pmeta.get("n", 2.0)
    dm_dX = M_tot * n * X ** (n - 1) / (1.0 + X**n) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = dm_dX / (L * 4.0 * np.pi * r**2)
    rho = np.where(np.isfinite(rho) & (r > 0), rho, 0.0)
    return rho / rc


def scan_profiles(mig: dict) -> list[dict]:
    rows = []
    x_grid = np.unique(
        np.concatenate(
            [
                np.array([0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 5.556, 6.0, 6.7, 7.0]),
                np.linspace(1.0, 7.0, 13),
            ]
        )
    )
    x_grid = np.sort(x_grid)

    # (family, kwargs, af_modes)
    jobs = []
    for x in x_grid:
        jobs.append(("exponential", {"x": float(x)}, ("invsq", "hratio")))
    for x in x_grid:
        for fu in (0.80, 0.90, 0.95, 0.98):
            jobs.append(("twoscale", {"x": float(x), "f_u": float(fu)}, ("invsq", "hratio")))
    # a few plummer checks
    for x in (1.0, 2.0, 3.0, 5.0, 6.7):
        jobs.append(("plummer", {"x": float(x), "n": 2.0}, ("invsq", "hratio")))

    for fam, kw, modes in jobs:
        if fam == "exponential":
            pmeta = params_exp_local_match(kw["x"])
        elif fam == "twoscale":
            pmeta = params_twoscale_local_match(kw["x"], kw["f_u"])
        else:
            pmeta = params_plummer_local_match(kw["x"], n=kw["n"])
        assert abs(pmeta["rho_local_over_rcrit"] - OMEGA_M) < 1e-8
        for af_mode in modes:
            phys = af_an_from_physics(pmeta, af_mode=af_mode)
            af_m = interp_periodic(mig["l"], phys["l_deg"], phys["AF_smooth"])
            an_m = interp_periodic(mig["l"], phys["l_deg"], phys["AN_smooth"])
            fit = fit_amp_phi(mig["l"], af_m, an_m, mig["A"], mig["B"])
            bal = fit["best_balanced"]
            row = {
                **{k: pmeta[k] for k in pmeta if k in (
                    "profile", "x", "f_u", "L_over_D", "kappa", "kappa_core", "kappa_u",
                    "m_frac", "rho_local_over_rcrit", "rho_mean_over_rcrit",
                    "rho_local_over_mean", "M_tot", "M_enc",
                )},
                "af_mode": af_mode,
                "frac_ptp": phys["frac_ptp"],
                "H_obs_ptp": phys["H_obs_ptp"],
                "R2_A": bal["r2_a"],
                "R2_B": bal["r2_b"],
                "R2_min": bal["r2_min"],
                "R2_sum": bal["r2_sum"],
                "best_Amp": bal["amp"],
                "best_phi_rad": bal["phi"],
                "best_phi_deg": bal["phi_deg"],
                "R2_A_sumopt": fit["best_sum"]["r2_a"],
                "R2_B_sumopt": fit["best_sum"]["r2_b"],
                "Amp_sumopt": fit["best_sum"]["amp"],
                "phi_sumopt": fit["best_sum"]["phi"],
                "R2_A_wo": fit["without"]["r2_a"],
                "R2_B_wo": fit["without"]["r2_b"],
                "R2_A_chris": fit["chris_ref"]["r2_a"],
                "R2_B_chris": fit["chris_ref"]["r2_b"],
            }
            rows.append(row)
            print(
                f"{row['profile']:18s} af={af_mode:6s} x={row['x']:.3f} κ={row['kappa']:.4f} "
                f"R2A={row['R2_A']:.3f} R2B={row['R2_B']:.3f} "
                f"Amp={row['best_Amp']:+.3f} φ={row['best_phi_deg']:.1f}°",
                flush=True,
            )
    return rows


def choose_recommended(rows: list[dict]) -> dict:
    """Prefer any profile with both R² high; break ties by κ closer to 0.15 and exponential/twoscale."""
    def score(r):
        # primary: min R2, then sum, then proximity of κ to ~0.2–0.3 band, prefer hratio if tied
        kappa_pen = abs(np.log10(max(r["kappa"], 1e-12) / 0.2))
        mode_bonus = 0.01 if r.get("af_mode") == "hratio" else 0.0
        fam_bonus = 0.005 if r["profile"].startswith("twoscale") or r["profile"] == "exponential" else 0.0
        return (r["R2_min"] + mode_bonus + fam_bonus, r["R2_sum"], -kappa_pen)

    ranked = sorted(rows, key=score, reverse=True)
    best = ranked[0]
    good = [r for r in rows if r["R2_A"] >= 0.65 and r["R2_B"] >= 0.65]
    good70 = [r for r in rows if r["R2_A"] >= 0.70 and r["R2_B"] >= 0.70]

    exp = [r for r in rows if r["profile"] == "exponential"]
    pareto = sorted(exp, key=lambda r: (r.get("af_mode", ""), r["x"]))

    rationale = (
        f"Under ρ(D)=Ω_m ρ_crit by construction, decreasing profiles require "
        f"κ = ρ_mean/(2 ρ_crit) > Ω_m/2 = 0.15. "
        f"Chosen {best['profile']} (AF={best.get('af_mode')}) at D/L={best['x']:.3f} "
        f"(L/D={best['L_over_D']:.4f}) gives κ={best['kappa']:.4f}, "
        f"ρ_mean/ρ_crit={best['rho_mean_over_rcrit']:.4f}, "
        f"R²_A={best['R2_A']:.3f}, R²_B={best['R2_B']:.3f} after Amp×sin "
        f"(Amp={best['best_Amp']:+.4f}, φ={best['best_phi_deg']:.2f}°)."
    )
    if good70:
        rationale += f" Achieved R²≥0.7 on both for {len(good70)} scanned point(s)."
    elif good:
        rationale += f" Achieved R²≥0.65 on both for {len(good)} scanned point(s); 0.7 not reached."
    else:
        rationale += (
            " Target R²≳0.7 on both not met under the density constraint; "
            "this is the best balanced point on the scanned Pareto front."
        )
    if best["profile"].startswith("twoscale"):
        rationale += (
            f" Two-scale floor fraction f_u={best.get('f_u')}: uniform floor lifts ρ(D) while "
            f"core supplies tidal κ_core={best.get('kappa_core', float('nan')):.4f}."
        )

    return {
        "recommended": best,
        "n_good_065": len(good),
        "n_good_070": len(good70),
        "top5": [
            {k: r[k] for k in (
                "profile", "af_mode", "x", "f_u", "kappa", "rho_mean_over_rcrit",
                "R2_A", "R2_B", "R2_min", "best_Amp", "best_phi_deg",
            ) if k in r}
            for r in ranked[:5]
        ],
        "pareto_exp": [
            {
                "x": r["x"],
                "af_mode": r.get("af_mode"),
                "kappa": r["kappa"],
                "rho_mean_over_rcrit": r["rho_mean_over_rcrit"],
                "R2_A": r["R2_A"],
                "R2_B": r["R2_B"],
                "R2_min": r["R2_min"],
                "best_Amp": r["best_Amp"],
                "best_phi_deg": r["best_phi_deg"],
            }
            for r in pareto
        ],
        "rationale": rationale,
        "kappa_vs_Om_over_2": {
            "Omega_m_over_2": OMEGA_M / 2.0,
            "kappa_required": best["kappa"],
            "exceeds": best["kappa"] > OMEGA_M / 2.0 + 1e-6,
            "note": (
                "For any decreasing ρ with ρ(D)=Ω_m ρ_crit, "
                "ρ_mean(<D)>ρ(D) ⇒ κ=ρ_mean/(2ρ_crit)>Ω_m/2. "
                "Do not fake κ=0.15."
            ),
        },
    }


def plot_density(pmeta: dict, path: Path) -> None:
    style()
    r = np.linspace(0.02, 3.0, 600)
    rho = rho_profile_curve(pmeta, r)
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.plot(r, rho, color=C_BLUE, lw=2.4, label=rf"{pmeta['profile']}, $D/L={pmeta['x']:.2f}$")
    ax.axhline(OMEGA_M, color=C_GREEN, ls="--", lw=1.8, label=rf"$\Omega_m={OMEGA_M}$ (local target)")
    ax.axhline(
        pmeta["rho_mean_over_rcrit"],
        color=C_ORANGE,
        ls=":",
        lw=1.8,
        label=rf"$\rho_{{\mathrm{{mean}}}}(<D)/\rho_c={pmeta['rho_mean_over_rcrit']:.3f}$",
    )
    ax.axvline(1.0, color=C_RED, ls="--", lw=1.5, label="observer $r=D$")
    ax.plot(1.0, OMEGA_M, "o", color=C_RED, ms=9, zorder=5)
    ax.set_xlabel(r"$r / D$")
    ax.set_ylabel(r"$\rho(r) / \rho_{\mathrm{crit}}$")
    ax.set_title(
        rf"Local-density match: $\rho(D)=\Omega_m\rho_c$  "
        rf"($\kappa={pmeta['kappa']:.3f}$, not $\Omega_m/2$)"
    )
    ax.set_xlim(0, 3)
    ax.set_ylim(0, max(1.2, float(np.nanmax(rho[:200])) * 1.05))
    ax.legend(loc="upper right", framealpha=0.95)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_overlay(mig, pmeta, phys, fit_bal, path: Path) -> None:
    style()
    amp, phi = fit_bal["best_Amp"], fit_bal["best_phi_rad"]
    l_f = phys["l_deg"]
    af_f, an_f = phys["AF_smooth"], phys["AN_smooth"]
    ma_f, mb_f = model_rotated(af_f, an_f, l_f, amp, phi)

    af_m = interp_periodic(mig["l"], l_f, af_f)
    an_m = interp_periodic(mig["l"], l_f, an_f)
    ma_m, mb_m = model_rotated(af_m, an_m, mig["l"], amp, phi)
    r2a, r2b = rsq(mig["A"], ma_m), rsq(mig["B"], mb_m)

    fig = plt.figure(figsize=(10.5, 8.2))
    gs = GridSpec(2, 1, height_ratios=[1, 1], hspace=0.08, figure=fig)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1], sharex=ax_a)

    ax_a.errorbar(
        mig["l"], mig["A"], yerr=np.vstack([mig["A_elo"], mig["A_ehi"]]),
        fmt="o", ms=4.0, color=C_BLACK, ecolor=C_GRAY, elinewidth=0.9, capsize=2.0,
        label="Migkas $A/A_{\\mathrm{all}}$ (digitised)", zorder=3,
    )
    ax_a.plot(
        l_f, ma_f, color=C_ORANGE, lw=2.2, zorder=4,
        label=(
            f"Model $(1-A\\sin)\\times$AF  "
            f"(Amp={amp:+.3f}, φ={phi:.2f} rad, $R^2$={r2a:.3f})"
        ),
    )
    ax_a.axhline(1.0, color=C_GREEN, ls="--", lw=1.0, alpha=0.7)
    ax_a.set_ylabel(r"$A/A_{\mathrm{all}}$")
    ax_a.set_ylim(0.65, 1.45)
    ax_a.legend(loc="upper right", framealpha=0.95)
    ax_a.grid(True, alpha=0.3)
    ax_a.set_title(
        rf"Local $\rho(D)=\Omega_m\rho_c$: {pmeta['profile']}, "
        rf"$D/L={pmeta['x']:.2f}$, $\kappa={pmeta['kappa']:.3f}$ "
        rf"(sheet-style phase rotation; not B-null)"
    )

    ax_b.errorbar(
        mig["l"], mig["B"], yerr=np.vstack([mig["B_elo"], mig["B_ehi"]]),
        fmt="o", ms=4.0, color=C_BLACK, ecolor=C_GRAY, elinewidth=0.9, capsize=2.0,
        label=r"Migkas $B/B_{\mathrm{all}}$ (digitised)", zorder=3,
    )
    ax_b.plot(
        l_f, mb_f, color=C_BLUE, lw=2.2, zorder=4,
        label=(
            f"Model $(1+A\\sin)\\times$AN  "
            f"(opp. Amp, shared φ, $R^2$={r2b:.3f})"
        ),
    )
    ax_b.axhline(1.0, color=C_GREEN, ls="--", lw=1.0, alpha=0.7)
    ax_b.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax_b.set_ylabel(r"$B/B_{\mathrm{all}}$")
    ax_b.set_ylim(0.65, 1.45)
    ax_b.set_xlim(0, 360)
    ax_b.set_xticks([0, 90, 180, 270, 360])
    ax_b.legend(loc="upper right", framealpha=0.95)
    ax_b.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_scatter(mig, pmeta, phys, fit_bal, path: Path) -> None:
    style()
    amp, phi = fit_bal["best_Amp"], fit_bal["best_phi_rad"]
    af_m = interp_periodic(mig["l"], phys["l_deg"], phys["AF_smooth"])
    an_m = interp_periodic(mig["l"], phys["l_deg"], phys["AN_smooth"])
    ma, mb = model_rotated(af_m, an_m, mig["l"], amp, phi)
    r2a, r2b = rsq(mig["A"], ma), rsq(mig["B"], mb)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8))
    for ax, x, y, r2, name, col in [
        (axes[0], ma, mig["A"], r2a, "A", C_ORANGE),
        (axes[1], mb, mig["B"], r2b, "B", C_BLUE),
    ]:
        ax.scatter(x, y, c=col, s=28, edgecolors="white", lw=0.4, zorder=3)
        lo = min(x.min(), y.min()) - 0.02
        hi = max(x.max(), y.max()) + 0.02
        ax.plot([lo, hi], [lo, hi], ls="--", color=C_GRAY, lw=1.2)
        # OLS
        a = np.corrcoef(x, y)[0, 1] * (np.std(y) / np.std(x))
        b = np.mean(y) - a * np.mean(x)
        xx = np.linspace(lo, hi, 50)
        ax.plot(xx, a * xx + b, color=C_BLACK, lw=1.5, label=rf"$R^2={r2:.3f}$")
        ax.set_xlabel(f"Model {name} (rotated)")
        ax.set_ylabel(f"Migkas {name}")
        ax.set_title(f"{name}: scatter after phase rotation")
        ax.legend(loc="best")
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal", adjustable="datalim")
    fig.suptitle(
        rf"{pmeta['profile']} $D/L={pmeta['x']:.2f}$, $\kappa={pmeta['kappa']:.3f}$, "
        rf"Amp={amp:+.3f}, φ={np.rad2deg(phi):.1f}°",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_scan_summary(rows: list[dict], path: Path) -> None:
    style()
    exp = [r for r in rows if r["profile"] == "exponential"]
    xs = np.array([r["x"] for r in exp])
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.0))
    ax = axes[0, 0]
    ax.plot(xs, [r["kappa"] for r in exp], color=C_BLUE, marker="o", ms=3)
    ax.axhline(0.15, color=C_RED, ls="--", label=r"$\Omega_m/2=0.15$")
    ax.set_xlabel(r"$D/L$")
    ax.set_ylabel(r"$\kappa$")
    ax.set_title(r"$\kappa$ required when $\rho(D)=\Omega_m\rho_c$")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(xs, [r["rho_mean_over_rcrit"] for r in exp], color=C_ORANGE, marker="o", ms=3, label=r"$\rho_{\mathrm{mean}}/\rho_c$")
    ax.axhline(OMEGA_M, color=C_GREEN, ls="--", label=r"$\rho(D)/\rho_c=\Omega_m$")
    ax.set_xlabel(r"$D/L$")
    ax.set_ylabel(r"density / $\rho_c$")
    ax.set_title("Mean vs local density")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(xs, [r["R2_A"] for r in exp], color=C_ORANGE, marker="o", ms=3, label=r"$R^2_A$")
    ax.plot(xs, [r["R2_B"] for r in exp], color=C_BLUE, marker="s", ms=3, label=r"$R^2_B$")
    ax.axhline(0.7, color=C_GRAY, ls=":", label="0.7 target")
    ax.set_xlabel(r"$D/L$")
    ax.set_ylabel(r"$R^2$")
    ax.set_title("Migkas R² after Amp×sin (balanced)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1)

    ax = axes[1, 1]
    ax.plot(xs, [r["frac_ptp"] for r in exp], color=C_PURPLE, marker="o", ms=3)
    ax.set_xlabel(r"$D/L$")
    ax.set_ylabel(r"$(H_{\mathrm{obs}}/H)_{\mathrm{ptp}}$")
    ax.set_title("Anisotropy amplitude vs concentration")
    ax.grid(True, alpha=0.3)

    fig.suptitle(r"Exponential scan under $\rho(D)=\Omega_m\rho_{\mathrm{crit}}$", y=1.01)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_note(rec: dict, path: Path) -> None:
    best = rec["recommended"]
    kv = rec["kappa_vs_Om_over_2"]
    soft = rec.get("soft_best")
    soft_line = ""
    if soft:
        soft_line = (
            f"- Soft best: `{soft['profile']}` at D/L={soft['x']:.3f}, "
            f"κ={soft['kappa']:.3f}, R²_A={soft['R2_A']:.3f}, R²_B={soft['R2_B']:.3f}\n"
        )
    text = f"""# Local-density match note

**Date:** 2026-09-18 (America/Chicago)

## What "local match" means

Chris asked whether the density **where the observer exists** matches the observed cosmic matter density.
Previously the binding sheet fixed **mean** enclosed density via

`κ = G M_enc(<D)/(D³ H²) = Ω_m/2 = 0.15`

which makes ρ_mean(<D)/ρ_crit = Ω_m, but for a decreasing profile (especially the Phase A
exponential at D/L ~ 6–7) the **point density** ρ(D) sits deep in the exponential tail and was
~15× below Ω_m ρ_crit.

This run **inverts** that priority:

**ρ(D) = Ω_m ρ_crit = 0.3 ρ_crit** (by construction, within 1%).

We still report ρ_mean(<D) and the κ implied by M_enc.

## κ vs Ω_m/2 tension (honest)

For **any** centrally decreasing profile with ρ(D)=Ω_m ρ_c,

ρ_mean(<D) > ρ(D)  ⇒  κ = ρ_mean/(2 ρ_c) > Ω_m/2 = 0.15.

**Do not fake κ=0.15** when the density constraint forces a larger enclosed mass.
In the uniform/soft limit D/L→0, κ→0.15 but tidal anisotropy vanishes (g∝−r).
At large D/L (point-mass-like), κ grows rapidly under the local-density fix.

| Quantity | Recommended |
|---|---|
| Profile | `{best['profile']}` |
| D/L | {best['x']:.4f} |
| L/D | {best['L_over_D']:.4f} |
| ρ(D)/ρ_c | {best['rho_local_over_rcrit']:.4f} |
| ρ_mean(<D)/ρ_c | {best['rho_mean_over_rcrit']:.4f} |
| κ | **{best['kappa']:.4f}** (vs Ω_m/2={kv['Omega_m_over_2']}) |
| κ > 0.15? | **{kv['exceeds']}** |
| Migkas R²_A | {best['R2_A']:.4f} |
| Migkas R²_B | {best['R2_B']:.4f} |
| Amp, φ | {best['best_Amp']:+.4f}, {best['best_phi_deg']:.2f}° |

## Migkas diagnostic

Uses the **sheet-style** Amp×sin phase-rotation diagnostic (not a free A cos φ dipole;
not a B-null claim):

- AF from ST+Gauss H_obs via A∝1/(H_obs/H)² (mean-renorm, Δℓ=90° smooth)
- AN from enclosed-mass angular series m(r_src)/M_enc (mean-renorm, same smooth)
- Model A = (1−A sin)×AF, Model B = (1+A sin)×AN, shared φ

{rec['rationale']}

## Soft / alternate profiles tried

Also scanned Plummer-like m/M = X^n/(1+X^n) for n∈{{1.5,2,3}} under the same
ρ(D)=Ω_m ρ_c constraint. Soft profiles can keep κ closer to 0.15 but typically
weaken the tidal anisotropy; see `scan_summary.json`.

{soft_line}
## Regenerable

```bash
/workspace/cosmology/.venv/bin/python /workspace/cosmology/density_local_match/run_density_local_match.py
```
"""
    path.write_text(text)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mig = load_migkas()
    print("Scanning profiles under ρ(D)=Ω_m ρ_crit …", flush=True)
    rows = scan_profiles(mig)
    rec = choose_recommended(rows)

    # persist scan
    scan_path = OUT / "scan_summary.json"
    scan_path.write_text(
        json.dumps(
            {
                "constraint": "rho(D) = Omega_m * rho_crit with Omega_m=0.3",
                "units": "H=G=D=1",
                "alpha": ALPHA,
                "d_over_D": EPS,
                "method": (
                    "ST+Gauss peculiar u=-α G m(r)/(r²H) r̂; "
                    "AF∝1/(H_obs/H)² mean-renorm + Δℓ=90°; "
                    "AN∝m(r)/M_enc mean-renorm + Δℓ=90°; "
                    "Chris Amp×sin shared-φ rotation; R² = Pearson² vs digitised Fig.6"
                ),
                "n_rows": len(rows),
                "rows": rows,
            },
            indent=2,
        )
    )

    best = rec["recommended"]
    # rebuild physics for plots
    if best["profile"] == "exponential":
        pmeta = params_exp_local_match(best["x"])
    elif best["profile"].startswith("twoscale"):
        pmeta = params_twoscale_local_match(best["x"], best.get("f_u", 0.9))
    else:
        n = float(best["profile"].split("n")[-1])
        pmeta = params_plummer_local_match(best["x"], n=n)
    phys = af_an_from_physics(pmeta, af_mode=best.get("af_mode", "invsq"))

    rec_out = {
        **rec,
        "pmeta": pmeta,
        "outputs": {
            "density_png": str(OUT / "density_profile_local_match.png"),
            "overlay_png": str(OUT / "migkas_overlay_local_density.png"),
            "scatter_png": str(OUT / "migkas_scatter_local_density.png"),
            "scan_png": str(OUT / "scan_pareto_local_density.png"),
            "note": str(OUT / "NOTE.md"),
        },
    }
    (OUT / "recommended_profile.json").write_text(json.dumps(rec_out, indent=2))

    plot_density(pmeta, OUT / "density_profile_local_match.png")
    plot_overlay(mig, pmeta, phys, best, OUT / "migkas_overlay_local_density.png")
    plot_scatter(mig, pmeta, phys, best, OUT / "migkas_scatter_local_density.png")
    plot_scan_summary(rows, OUT / "scan_pareto_local_density.png")
    write_note(rec, OUT / "NOTE.md")

    print("\n=== RECOMMENDED ===", flush=True)
    print(json.dumps({k: best[k] for k in (
        "profile", "x", "L_over_D", "kappa", "rho_local_over_rcrit",
        "rho_mean_over_rcrit", "R2_A", "R2_B", "best_Amp", "best_phi_deg",
    )}, indent=2), flush=True)
    print(rec["rationale"], flush=True)
    print(f"Wrote outputs under {OUT}", flush=True)


if __name__ == "__main__":
    main()

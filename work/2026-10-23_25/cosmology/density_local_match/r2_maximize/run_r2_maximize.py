#!/usr/bin/env python3
"""
Aggressive R² maximisation under ρ(D)=Ω_m ρ_crit (Ω_m=0.3).

Extends density_local_match: scan profile family, D/L, f_u, α, d/D, AF mode,
Amp/φ — including independent Amp_A/Amp_B and separate φ for A-only / B-only.

Reports three optima:
  1) Best A alone (max R²_A)
  2) Best B alone (max R²_B)
  3) Best pair: max min(R²_A,R²_B) and max mean/sum

Do NOT fake κ=0.15 — report actual κ from M_enc.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
ROOT = Path("/workspace/cosmology")
OUT = HERE  # write into r2_maximize/
CSV = ROOT / "migkas" / "fig6_digitized.csv"

# Load sibling module without requiring package install
_spec = importlib.util.spec_from_file_location(
    "run_density_local_match", PARENT / "run_density_local_match.py"
)
_mod = importlib.util.module_from_spec(_spec)
sys.modules["run_density_local_match"] = _mod
_spec.loader.exec_module(_mod)

OMEGA_M = _mod.OMEGA_M
OMEGA = _mod.OMEGA
H = _mod.H
G = _mod.G
D = _mod.D
C_BLUE = _mod.C_BLUE
C_ORANGE = _mod.C_ORANGE
C_GREEN = _mod.C_GREEN
C_BLACK = _mod.C_BLACK
C_GRAY = _mod.C_GRAY
C_RED = _mod.C_RED
C_PURPLE = _mod.C_PURPLE

params_exp_local_match = _mod.params_exp_local_match
params_twoscale_local_match = _mod.params_twoscale_local_match
params_plummer_local_match = _mod.params_plummer_local_match
m_enc_at_r = _mod.m_enc_at_r
boxcar_dl90 = _mod.boxcar_dl90
rsq = _mod.rsq
sine_mod = _mod.sine_mod
interp_periodic = _mod.interp_periodic
load_migkas = _mod.load_migkas
style = _mod.style
rho_crit = _mod.rho_crit


def peculiar_u(r_vec, pmeta, alpha: float):
    r = np.linalg.norm(r_vec, axis=-1)
    r_safe = np.maximum(r, 1e-30)
    m = m_enc_at_r(r_safe, pmeta)
    pref = -alpha * G * m / (r_safe**3 * H)
    return pref[..., None] * r_vec


def H_obs_circle(pmeta, eps: float, alpha: float, n_phi: int = 721):
    l_deg = np.linspace(0.0, 360.0, n_phi, endpoint=False)
    phi = np.deg2rad(l_deg)
    n_hat = np.stack([np.cos(phi), np.sin(phi), np.zeros_like(phi)], axis=-1)
    d = eps * D
    D_vec = np.array([D, 0.0, 0.0])
    r_vec = D_vec + d * n_hat
    u_r = peculiar_u(r_vec, pmeta, alpha=alpha)
    u_D = peculiar_u(D_vec.reshape(1, 3), pmeta, alpha=alpha)[0]
    delta_v = np.sum((u_r - u_D) * n_hat, axis=-1)
    H_obs = H + delta_v / d
    m_src = m_enc_at_r(np.linalg.norm(r_vec, axis=-1), pmeta)
    return l_deg, H_obs, m_src


def af_an_from_physics(pmeta, eps: float, alpha: float, af_mode: str = "hratio"):
    l_deg, H_obs, m_src = H_obs_circle(pmeta, eps=eps, alpha=alpha)
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


def model_AB(af, an, l_deg, amp_a, amp_b, phi):
    """Independent Amp_A, Amp_B with shared φ.
    A = (1 - Amp_A sin)×AF;  B = (1 - Amp_B sin)×AN.
    Shared-|Amp| opposite-sign case is Amp_B = -Amp_A.
    """
    s = sine_mod(l_deg, 1.0, phi)
    return (1.0 - amp_a * s) * af, (1.0 - amp_b * s) * an


def model_A_only(af, l_deg, amp_a, phi_a):
    return (1.0 - sine_mod(l_deg, amp_a, phi_a)) * af


def model_B_only(an, l_deg, amp_b, phi_b):
    return (1.0 - sine_mod(l_deg, amp_b, phi_b)) * an


def _amp_grid(lo=0.02, hi=0.20, n=19, include_neg=True, include_zero=True):
    pos = np.linspace(lo, hi, n)
    parts = []
    if include_zero:
        parts.append(np.array([0.0]))
    if include_neg:
        parts.append(-pos)
    parts.append(pos)
    return np.unique(np.concatenate(parts))


def _pearson_r2_batch(y, X):
    """R² of y vs each column of X (n × k). y shape (n,)."""
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    y0 = y - y.mean()
    X0 = X - X.mean(axis=0, keepdims=True)
    den_y = np.sqrt(np.sum(y0 * y0))
    den_x = np.sqrt(np.sum(X0 * X0, axis=0))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (X0.T @ y0) / (den_x * den_y + 1e-30)
        r2 = r * r
    r2 = np.where(np.isfinite(r2), r2, 0.0)
    return r2


def fit_shared_opp(l, af, an, a_obs, b_obs, amp_grid=None, n_phi=181):
    """Shared |Amp| opposite sign (Amp_B=-Amp_A), shared φ. Vectorized."""
    if amp_grid is None:
        amp_grid = _amp_grid(0.02, 0.20, 37)
    amps = np.asarray(amp_grid, dtype=float)
    phis = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    # s[phi, ell]
    s = np.sin(OMEGA * np.asarray(l, dtype=float)[None, :] + phis[:, None])
    af = np.asarray(af, dtype=float)
    an = np.asarray(an, dtype=float)
    # For each (phi, amp): ma = af*(1 - amp*s) → shape (n_phi, n_amp, n_ell) too big
    # Process per phi (n_phi ~ 120–360, n_amp ~ 40, n_ell ~ 70) → OK if we loop phi lightly
    best_min = {"r2_min": -1.0, "r2_sum": -1.0}
    best_mean = {"r2_sum": -1.0, "r2_min": -1.0}
    for i, phi in enumerate(phis):
        si = s[i]
        # ma[:, ell] for each amp: af - amp*(af*si)
        base_a = af * si
        base_b = an * si
        Ma = af[None, :] - amps[:, None] * base_a[None, :]
        Mb = an[None, :] - (-amps)[:, None] * base_b[None, :]  # Amp_B = -Amp_A
        # Wait: B = (1 - Amp_B * s)*AN with Amp_B=-amp ⇒ (1 - (-amp)*s)*AN = (1+amp*s)*AN
        # = AN - (-amp)*(AN*s) = AN + amp*(AN*s)
        Mb = an[None, :] - (-amps)[:, None] * base_b[None, :]
        ra = _pearson_r2_batch(a_obs, Ma.T)
        rb = _pearson_r2_batch(b_obs, Mb.T)
        ssum = ra + rb
        mn = np.minimum(ra, rb)
        j_min = int(np.argmax(mn + 1e-9 * ssum))
        j_mean = int(np.argmax(ssum + 1e-9 * mn))
        for j, which in ((j_min, "min"), (j_mean, "mean")):
            rec = {
                "mode": "shared_opp",
                "amp_a": float(amps[j]),
                "amp_b": float(-amps[j]),
                "phi": float(phi),
                "phi_deg": float(np.rad2deg(phi) % 360.0),
                "r2_a": float(ra[j]),
                "r2_b": float(rb[j]),
                "r2_sum": float(ssum[j]),
                "r2_mean": float(ssum[j] / 2.0),
                "r2_min": float(mn[j]),
            }
            if which == "min":
                if rec["r2_min"] > best_min["r2_min"] or (
                    abs(rec["r2_min"] - best_min["r2_min"]) < 1e-14
                    and rec["r2_sum"] > best_min["r2_sum"]
                ):
                    best_min = rec
            else:
                if rec["r2_sum"] > best_mean["r2_sum"] or (
                    abs(rec["r2_sum"] - best_mean["r2_sum"]) < 1e-14
                    and rec["r2_min"] > best_mean["r2_min"]
                ):
                    best_mean = rec
    return {"best_min": best_min, "best_mean": best_mean}


def fit_indep_amps_shared_phi(l, af, an, a_obs, b_obs, amp_grid=None, n_phi=181):
    """Independent Amp_A, Amp_B; shared φ. Vectorized per φ."""
    if amp_grid is None:
        amp_grid = _amp_grid(0.02, 0.20, 37)
    amps = np.asarray(amp_grid, dtype=float)
    phis = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    af = np.asarray(af, dtype=float)
    an = np.asarray(an, dtype=float)
    best_min = {"r2_min": -1.0, "r2_sum": -1.0}
    best_mean = {"r2_sum": -1.0, "r2_min": -1.0}
    for phi in phis:
        si = np.sin(OMEGA * np.asarray(l, dtype=float) + float(phi))
        Ma = af[None, :] - amps[:, None] * (af * si)[None, :]
        Mb = an[None, :] - amps[:, None] * (an * si)[None, :]
        ra = _pearson_r2_batch(a_obs, Ma.T)
        rb = _pearson_r2_batch(b_obs, Mb.T)
        ja = int(np.argmax(ra))
        jb = int(np.argmax(rb))
        best_ra, best_aa = float(ra[ja]), float(amps[ja])
        best_rb, best_ab = float(rb[jb]), float(amps[jb])
        ssum = best_ra + best_rb
        mn = min(best_ra, best_rb)
        rec = {
            "mode": "indep_amps_shared_phi",
            "amp_a": best_aa,
            "amp_b": best_ab,
            "phi": float(phi),
            "phi_deg": float(np.rad2deg(phi) % 360.0),
            "r2_a": best_ra,
            "r2_b": best_rb,
            "r2_sum": float(ssum),
            "r2_mean": float(ssum / 2.0),
            "r2_min": float(mn),
        }
        if mn > best_min["r2_min"] or (
            abs(mn - best_min["r2_min"]) < 1e-14 and ssum > best_min["r2_sum"]
        ):
            best_min = rec
        if ssum > best_mean["r2_sum"] or (
            abs(ssum - best_mean["r2_sum"]) < 1e-14 and mn > best_mean["r2_min"]
        ):
            best_mean = rec
    return {"best_min": best_min, "best_mean": best_mean}


def fit_solo(l, base, obs, amp_grid=None, n_phi=181, which="A"):
    """Independent Amp and φ for one series. Vectorized."""
    if amp_grid is None:
        amp_grid = _amp_grid(0.02, 0.25, 47)
    amps = np.asarray(amp_grid, dtype=float)
    phis = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    base = np.asarray(base, dtype=float)
    best = {"r2": -1.0}
    for phi in phis:
        si = np.sin(OMEGA * np.asarray(l, dtype=float) + float(phi))
        M = base[None, :] - amps[:, None] * (base * si)[None, :]
        r2 = _pearson_r2_batch(obs, M.T)
        j = int(np.argmax(r2))
        if float(r2[j]) > best["r2"]:
            best = {
                "which": which,
                "amp": float(amps[j]),
                "phi": float(phi),
                "phi_deg": float(np.rad2deg(phi) % 360.0),
                "r2": float(r2[j]),
            }
    return best



def make_pmeta(fam: str, x: float, f_u: float | None = None, n: float = 2.0):
    if fam == "exponential":
        return params_exp_local_match(x)
    if fam == "twoscale":
        return params_twoscale_local_match(x, f_u)
    if fam == "plummer":
        return params_plummer_local_match(x, n=n)
    raise ValueError(fam)


def physics_row(pmeta, af_mode, alpha, eps, mig, refine_phi=False):
    """Evaluate one physics point; return pair + solo fits."""
    assert abs(pmeta["rho_local_over_rcrit"] - OMEGA_M) < 0.01 * OMEGA_M + 1e-9
    phys = af_an_from_physics(pmeta, eps=eps, alpha=alpha, af_mode=af_mode)
    af_m = interp_periodic(mig["l"], phys["l_deg"], phys["AF_smooth"])
    an_m = interp_periodic(mig["l"], phys["l_deg"], phys["AN_smooth"])

    n_phi = 361 if refine_phi else 121
    amp_g = _amp_grid(0.02, 0.22, 41 if refine_phi else 21)

    shared = fit_shared_opp(mig["l"], af_m, an_m, mig["A"], mig["B"], amp_grid=amp_g, n_phi=n_phi)
    indep = fit_indep_amps_shared_phi(
        mig["l"], af_m, an_m, mig["A"], mig["B"], amp_grid=amp_g, n_phi=n_phi
    )
    # Pick better of shared-opp vs indep-amps for each criterion
    cand_min = [shared["best_min"], indep["best_min"]]
    cand_mean = [shared["best_mean"], indep["best_mean"]]
    pair_min = max(cand_min, key=lambda r: (r["r2_min"], r["r2_sum"]))
    pair_mean = max(cand_mean, key=lambda r: (r["r2_sum"], r["r2_min"]))

    solo_a = fit_solo(mig["l"], af_m, mig["A"], amp_grid=amp_g, n_phi=n_phi, which="A")
    solo_b = fit_solo(mig["l"], an_m, mig["B"], amp_grid=amp_g, n_phi=n_phi, which="B")

    # Also report shared-opp R² at solo-A φ (diagnostic) — skip for speed

    base = {
        "profile": pmeta["profile"],
        "fam": pmeta["profile"].split("_")[0] if "_" in pmeta["profile"] else pmeta["profile"],
        "x": float(pmeta["x"]),
        "f_u": float(pmeta["f_u"]) if "f_u" in pmeta else None,
        "L_over_D": float(pmeta["L_over_D"]),
        "kappa": float(pmeta["kappa"]),
        "kappa_core": float(pmeta["kappa_core"]) if "kappa_core" in pmeta else None,
        "kappa_u": float(pmeta["kappa_u"]) if "kappa_u" in pmeta else None,
        "rho_local_over_rcrit": float(pmeta["rho_local_over_rcrit"]),
        "rho_mean_over_rcrit": float(pmeta["rho_mean_over_rcrit"]),
        "M_enc": float(pmeta["M_enc"]),
        "af_mode": af_mode,
        "alpha": float(alpha),
        "eps": float(eps),
        "frac_ptp": phys["frac_ptp"],
        "H_obs_ptp": phys["H_obs_ptp"],
        # pair (min criterion winner)
        "pair_min_mode": pair_min["mode"],
        "R2_A_mincrit": pair_min["r2_a"],
        "R2_B_mincrit": pair_min["r2_b"],
        "R2_min": pair_min["r2_min"],
        "R2_mean_at_mincrit": pair_min["r2_mean"],
        "Amp_A_mincrit": pair_min["amp_a"],
        "Amp_B_mincrit": pair_min["amp_b"],
        "phi_mincrit_deg": pair_min["phi_deg"],
        "phi_mincrit_rad": pair_min["phi"],
        # pair (mean/sum criterion)
        "pair_mean_mode": pair_mean["mode"],
        "R2_A_meancrit": pair_mean["r2_a"],
        "R2_B_meancrit": pair_mean["r2_b"],
        "R2_sum": pair_mean["r2_sum"],
        "R2_mean": pair_mean["r2_mean"],
        "R2_min_at_meancrit": pair_mean["r2_min"],
        "Amp_A_meancrit": pair_mean["amp_a"],
        "Amp_B_meancrit": pair_mean["amp_b"],
        "phi_meancrit_deg": pair_mean["phi_deg"],
        "phi_meancrit_rad": pair_mean["phi"],
        # solo
        "R2_A_solo": solo_a["r2"],
        "Amp_A_solo": solo_a["amp"],
        "phi_A_solo_deg": solo_a["phi_deg"],
        "phi_A_solo_rad": solo_a["phi"],
        "R2_B_solo": solo_b["r2"],
        "Amp_B_solo": solo_b["amp"],
        "phi_B_solo_deg": solo_b["phi_deg"],
        "phi_B_solo_rad": solo_b["phi"],
        # convenience aliases matching baseline
        "R2_A": pair_min["r2_a"],
        "R2_B": pair_min["r2_b"],
        "best_Amp": pair_min["amp_a"],
        "best_phi_deg": pair_min["phi_deg"],
        "best_phi_rad": pair_min["phi"],
    }
    return base, phys, af_m, an_m


def coarse_jobs():
    jobs = []
    x_exp = np.unique(np.concatenate([
        np.array([0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 5.5, 6.0, 6.7, 7.0, 8.0]),
        np.linspace(1.0, 3.0, 9),
    ]))
    x_ts = np.unique(np.concatenate([
        np.linspace(3.0, 8.0, 21),
        np.array([5.2, 5.4, 5.6, 5.8, 5.9, 6.0, 6.2, 6.4, 6.6]),
    ]))
    f_us = np.array([0.70, 0.80, 0.85, 0.88, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98])
    alphas = np.array([0.5, 0.75, 1.0, 1.25, 1.5, 2.0])
    epss = np.array([0.05, 0.08, 0.10, 0.12, 0.15, 0.20])
    af_modes = ("hratio", "invsq")

    # Focus: baseline neighbourhood denser; global coarser on alpha/eps first
    # Stage strategy handled in main: here list all, main may subsample

    for x in x_exp:
        for af in af_modes:
            for a in alphas:
                for e in epss:
                    jobs.append(("exponential", float(x), None, af, float(a), float(e)))
    for x in x_ts:
        for fu in f_us:
            for af in af_modes:
                for a in alphas:
                    for e in epss:
                        jobs.append(("twoscale", float(x), float(fu), af, float(a), float(e)))
    # a few plummer
    for x in (1.0, 2.0, 3.0, 5.0):
        for af in af_modes:
            for a in (1.0, 1.5):
                for e in (0.1, 0.15):
                    jobs.append(("plummer", float(x), None, af, float(a), float(e)))
    return jobs


def stage1_jobs():
    """Coarser first pass — fewer alpha/eps combos."""
    jobs = []
    x_exp = np.array([0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0])
    x_ts = np.array([4.5, 5.2, 5.6, 5.9, 6.2, 6.6, 7.0, 7.5])
    f_us = np.array([0.85, 0.90, 0.92, 0.94, 0.95, 0.96, 0.97, 0.98])
    alphas = np.array([0.75, 1.0, 1.25, 1.5])
    epss = np.array([0.08, 0.10, 0.12, 0.15])
    # Exponential full
    for x in x_exp:
        for af in ("hratio", "invsq"):
            for a in alphas:
                for e in epss:
                    jobs.append(("exponential", float(x), None, af, float(a), float(e)))
    # Two-scale hratio full α×eps; invsq only near α=1,eps=0.1
    for x in x_ts:
        for fu in f_us:
            for a in alphas:
                for e in epss:
                    jobs.append(("twoscale", float(x), float(fu), "hratio", float(a), float(e)))
            for a in (1.0, 1.25):
                for e in (0.10, 0.12):
                    jobs.append(("twoscale", float(x), float(fu), "invsq", float(a), float(e)))
    for x in (1.5, 2.5, 5.0):
        for af in ("hratio", "invsq"):
            jobs.append(("plummer", float(x), None, af, 1.0, 0.10))
    return jobs



def refine_jobs_around(seeds: list[dict], n_each=5):
    """Tight local refine grids around top seeds (kept small on purpose)."""
    jobs = []
    seen = set()
    for s in seeds[:n_each]:
        x0 = s["x"]
        fu0 = s.get("f_u")
        a0 = s["alpha"]
        e0 = s["eps"]
        af = s["af_mode"]
        fam = "twoscale" if (fu0 is not None or str(s["profile"]).startswith("twoscale")) else (
            "plummer" if str(s["profile"]).startswith("plummer") else "exponential"
        )
        xs = np.unique(np.clip(np.linspace(x0 - 0.4, x0 + 0.4, 5), 0.5, 10.0))
        alphas = np.unique(np.clip(np.linspace(a0 - 0.25, a0 + 0.25, 5), 0.3, 3.0))
        epss = np.unique(np.clip(np.linspace(e0 - 0.03, e0 + 0.03, 3), 0.03, 0.30))
        if fam == "twoscale":
            fus = np.unique(np.clip(np.linspace(fu0 - 0.03, fu0 + 0.03, 5), 0.5, 0.995))
            for x in xs:
                for fu in fus:
                    for a in alphas:
                        for e in epss:
                            key = (fam, round(x, 4), round(fu, 4), af, round(a, 4), round(e, 4))
                            if key not in seen:
                                seen.add(key)
                                jobs.append((fam, float(x), float(fu), af, float(a), float(e)))
        else:
            for x in xs:
                for a in alphas:
                    for e in epss:
                        key = (fam, round(x, 4), None, af, round(a, 4), round(e, 4))
                        if key not in seen:
                            seen.add(key)
                            jobs.append((fam, float(x), None, af, float(a), float(e)))
    return jobs


def run_jobs(jobs, mig, refine_phi=False, label=""):
    rows = []
    t0 = time.time()
    n = len(jobs)
    for i, (fam, x, fu, af, alpha, eps) in enumerate(jobs):
        try:
            if fam == "twoscale":
                pmeta = make_pmeta(fam, x, f_u=fu)
            elif fam == "plummer":
                pmeta = make_pmeta(fam, x, n=2.0)
            else:
                pmeta = make_pmeta(fam, x)
            row, _, _, _ = physics_row(pmeta, af, alpha, eps, mig, refine_phi=refine_phi)
            rows.append(row)
        except Exception as ex:
            print(f"  skip {fam} x={x} fu={fu}: {ex}", flush=True)
            continue
        if (i + 1) % 25 == 0 or i + 1 == n:
            elapsed = time.time() - t0
            best_min = max((r["R2_min"] for r in rows), default=float("nan"))
            best_a = max((r["R2_A_solo"] for r in rows), default=float("nan"))
            best_b = max((r["R2_B_solo"] for r in rows), default=float("nan"))
            print(
                f"  [{label}] {i+1}/{n}  t={elapsed:.0f}s  "
                f"best_min={best_min:.4f}  best_Asolo={best_a:.4f}  best_Bsolo={best_b:.4f}",
                flush=True,
            )
    return rows


def plot_winner(mig, pmeta, phys, amp_a, amp_b, phi, r2a, r2b, title, path_overlay, path_scatter=None):
    style()
    l_f = phys["l_deg"]
    af_f, an_f = phys["AF_smooth"], phys["AN_smooth"]
    ma_f, mb_f = model_AB(af_f, an_f, l_f, amp_a, amp_b, phi)

    af_m = interp_periodic(mig["l"], l_f, af_f)
    an_m = interp_periodic(mig["l"], l_f, an_f)
    ma_m, mb_m = model_AB(af_m, an_m, mig["l"], amp_a, amp_b, phi)
    # recompute r2 for consistency
    r2a = rsq(mig["A"], ma_m)
    r2b = rsq(mig["B"], mb_m)

    fig = plt.figure(figsize=(10.5, 8.2))
    gs = GridSpec(2, 1, height_ratios=[1, 1], hspace=0.08, figure=fig)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1], sharex=ax_a)

    ax_a.errorbar(
        mig["l"], mig["A"], yerr=np.vstack([mig["A_elo"], mig["A_ehi"]]),
        fmt="o", ms=4.0, color=C_BLACK, ecolor=C_GRAY, elinewidth=0.9, capsize=2.0,
        label="Migkas $A/A_{\\mathrm{all}}$", zorder=3,
    )
    ax_a.plot(
        l_f, ma_f, color=C_ORANGE, lw=2.2, zorder=4,
        label=f"Model $(1-A_A\\sin)\\times$AF  (Amp_A={amp_a:+.3f}, φ={np.rad2deg(phi)%360:.1f}°, $R^2$={r2a:.3f})",
    )
    ax_a.axhline(1.0, color=C_GREEN, ls="--", lw=1.0, alpha=0.7)
    ax_a.set_ylabel(r"$A/A_{\mathrm{all}}$")
    ax_a.set_ylim(0.65, 1.45)
    ax_a.legend(loc="upper right", framealpha=0.95, fontsize=8)
    ax_a.grid(True, alpha=0.3)
    ax_a.set_title(title)

    ax_b.errorbar(
        mig["l"], mig["B"], yerr=np.vstack([mig["B_elo"], mig["B_ehi"]]),
        fmt="o", ms=4.0, color=C_BLACK, ecolor=C_GRAY, elinewidth=0.9, capsize=2.0,
        label=r"Migkas $B/B_{\mathrm{all}}$", zorder=3,
    )
    ax_b.plot(
        l_f, mb_f, color=C_BLUE, lw=2.2, zorder=4,
        label=f"Model $(1-A_B\\sin)\\times$AN  (Amp_B={amp_b:+.3f}, $R^2$={r2b:.3f})",
    )
    ax_b.axhline(1.0, color=C_GREEN, ls="--", lw=1.0, alpha=0.7)
    ax_b.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax_b.set_ylabel(r"$B/B_{\mathrm{all}}$")
    ax_b.set_ylim(0.65, 1.45)
    ax_b.set_xlim(0, 360)
    ax_b.set_xticks([0, 90, 180, 270, 360])
    ax_b.legend(loc="upper right", framealpha=0.95, fontsize=8)
    ax_b.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path_overlay)
    plt.close(fig)

    if path_scatter is not None:
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8))
        for ax, x, y, r2, name, col in [
            (axes[0], ma_m, mig["A"], r2a, "A", C_ORANGE),
            (axes[1], mb_m, mig["B"], r2b, "B", C_BLUE),
        ]:
            ax.scatter(x, y, c=col, s=28, edgecolors="white", lw=0.4, zorder=3)
            lo = min(x.min(), y.min()) - 0.02
            hi = max(x.max(), y.max()) + 0.02
            ax.plot([lo, hi], [lo, hi], ls="--", color=C_GRAY, lw=1.2)
            a = np.corrcoef(x, y)[0, 1] * (np.std(y) / np.std(x))
            b = np.mean(y) - a * np.mean(x)
            xx = np.linspace(lo, hi, 50)
            ax.plot(xx, a * xx + b, color=C_BLACK, lw=1.5, label=rf"$R^2={r2:.3f}$")
            ax.set_xlabel(f"Model {name}")
            ax.set_ylabel(f"Migkas {name}")
            ax.legend(loc="best")
            ax.grid(True, alpha=0.3)
            ax.set_aspect("equal", adjustable="datalim")
        fig.suptitle(title, y=1.02, fontsize=11)
        fig.tight_layout()
        fig.savefig(path_scatter, bbox_inches="tight")
        plt.close(fig)
    return r2a, r2b


def plot_solo_combined(mig, pmeta, phys, amp_a, phi_a, amp_b, phi_b, r2a, r2b, title, path):
    """Overlay with independent φ_A, φ_B (A-only / B-only winners shown together)."""
    style()
    l_f = phys["l_deg"]
    af_f, an_f = phys["AF_smooth"], phys["AN_smooth"]
    ma_f = model_A_only(af_f, l_f, amp_a, phi_a)
    mb_f = model_B_only(an_f, l_f, amp_b, phi_b)

    fig = plt.figure(figsize=(10.5, 8.2))
    gs = GridSpec(2, 1, height_ratios=[1, 1], hspace=0.08, figure=fig)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1], sharex=ax_a)
    ax_a.errorbar(
        mig["l"], mig["A"], yerr=np.vstack([mig["A_elo"], mig["A_ehi"]]),
        fmt="o", ms=4.0, color=C_BLACK, ecolor=C_GRAY, elinewidth=0.9, capsize=2.0,
        label="Migkas A", zorder=3,
    )
    ax_a.plot(
        l_f, ma_f, color=C_ORANGE, lw=2.2,
        label=f"A-only Amp={amp_a:+.3f}, φ={np.rad2deg(phi_a)%360:.1f}°, $R^2$={r2a:.3f}",
    )
    ax_a.axhline(1.0, color=C_GREEN, ls="--", lw=1.0, alpha=0.7)
    ax_a.set_ylabel(r"$A/A_{\mathrm{all}}$")
    ax_a.set_ylim(0.65, 1.45)
    ax_a.legend(loc="upper right", fontsize=8)
    ax_a.grid(True, alpha=0.3)
    ax_a.set_title(title)

    ax_b.errorbar(
        mig["l"], mig["B"], yerr=np.vstack([mig["B_elo"], mig["B_ehi"]]),
        fmt="o", ms=4.0, color=C_BLACK, ecolor=C_GRAY, elinewidth=0.9, capsize=2.0,
        label="Migkas B", zorder=3,
    )
    ax_b.plot(
        l_f, mb_f, color=C_BLUE, lw=2.2,
        label=f"B-only Amp={amp_b:+.3f}, φ={np.rad2deg(phi_b)%360:.1f}°, $R^2$={r2b:.3f}",
    )
    ax_b.axhline(1.0, color=C_GREEN, ls="--", lw=1.0, alpha=0.7)
    ax_b.set_xlabel(r"Galactic longitude $\ell$ (deg)")
    ax_b.set_ylabel(r"$B/B_{\mathrm{all}}$")
    ax_b.set_ylim(0.65, 1.45)
    ax_b.set_xlim(0, 360)
    ax_b.set_xticks([0, 90, 180, 270, 360])
    ax_b.legend(loc="upper right", fontsize=8)
    ax_b.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def rebuild_phys(row):
    fam = "twoscale" if row.get("f_u") is not None or str(row["profile"]).startswith("twoscale") else (
        "plummer" if str(row["profile"]).startswith("plummer") else "exponential"
    )
    if fam == "twoscale":
        pmeta = make_pmeta(fam, row["x"], f_u=row["f_u"])
    elif fam == "plummer":
        pmeta = make_pmeta(fam, row["x"], n=2.0)
    else:
        pmeta = make_pmeta(fam, row["x"])
    phys = af_an_from_physics(pmeta, eps=row["eps"], alpha=row["alpha"], af_mode=row["af_mode"])
    return pmeta, phys


def slim(row, keys=None):
    if keys is None:
        keys = [
            "profile", "af_mode", "x", "f_u", "alpha", "eps", "kappa",
            "rho_local_over_rcrit", "rho_mean_over_rcrit",
            "R2_min", "R2_A_mincrit", "R2_B_mincrit", "R2_mean", "R2_sum",
            "R2_A_meancrit", "R2_B_meancrit",
            "pair_min_mode", "Amp_A_mincrit", "Amp_B_mincrit", "phi_mincrit_deg",
            "pair_mean_mode", "Amp_A_meancrit", "Amp_B_meancrit", "phi_meancrit_deg",
            "R2_A_solo", "Amp_A_solo", "phi_A_solo_deg",
            "R2_B_solo", "Amp_B_solo", "phi_B_solo_deg",
            "frac_ptp",
        ]
    return {k: row[k] for k in keys if k in row and row[k] is not None}


def write_note(summary, path: Path):
    wA = summary["best_A_only"]
    wB = summary["best_B_only"]
    wPmin = summary["best_pair_min"]
    wPmean = summary["best_pair_mean"]
    text = f"""# R² maximise note (local-density constraint)

**Date:** 2026-09-18 (America/Chicago)

## Constraint (unchanged)

ρ(D) = Ω_m ρ_crit = **0.3 ρ_crit** (verified within 1%).
κ is **derived** from M_enc — **not** forced to 0.15.
For decreasing ρ with this local match, κ > Ω_m/2 = 0.15 always.

## Baseline to beat

Prior balanced (shared |Amp| opp. sign): **R²_A ≈ R²_B ≈ 0.685**, κ≈0.231, D/L≈5.9, f_u≈0.94, AF=`hratio`.

Chris sheet ballpark ~**0.74 / 0.73** (different normalisation / freer phase bookkeeping; not under this exact local-density fix).

## Three optima (this scan)

### 1) Best A alone (independent Amp_A, φ_A)

| | |
|---|---|
| R²_A | **{wA['R2_A_solo']:.4f}** |
| Profile | `{wA['profile']}` |
| D/L, f_u | {wA['x']:.4f}, {wA.get('f_u')} |
| α, d/D | {wA['alpha']:.3f}, {wA['eps']:.3f} |
| AF | `{wA['af_mode']}` |
| κ | **{wA['kappa']:.4f}** |
| ρ_mean/ρ_c | {wA['rho_mean_over_rcrit']:.4f} |
| Amp_A, φ_A | {wA['Amp_A_solo']:+.4f}, {wA['phi_A_solo_deg']:.2f}° |

### 2) Best B alone (independent Amp_B, φ_B)

| | |
|---|---|
| R²_B | **{wB['R2_B_solo']:.4f}** |
| Profile | `{wB['profile']}` |
| D/L, f_u | {wB['x']:.4f}, {wB.get('f_u')} |
| α, d/D | {wB['alpha']:.3f}, {wB['eps']:.3f} |
| AF | `{wB['af_mode']}` |
| κ | **{wB['kappa']:.4f}** |
| Amp_B, φ_B | {wB['Amp_B_solo']:+.4f}, {wB['phi_B_solo_deg']:.2f}° |

### 3a) Best pair — max min(R²_A, R²_B)

| | |
|---|---|
| R²_A / R²_B | **{wPmin['R2_A_mincrit']:.4f} / {wPmin['R2_B_mincrit']:.4f}** |
| R²_min | **{wPmin['R2_min']:.4f}** |
| Mode | `{wPmin['pair_min_mode']}` |
| Profile | `{wPmin['profile']}` |
| D/L, f_u | {wPmin['x']:.4f}, {wPmin.get('f_u')} |
| α, d/D | {wPmin['alpha']:.3f}, {wPmin['eps']:.3f} |
| AF | `{wPmin['af_mode']}` |
| κ | **{wPmin['kappa']:.4f}** |
| Amp_A, Amp_B, φ | {wPmin['Amp_A_mincrit']:+.4f}, {wPmin['Amp_B_mincrit']:+.4f}, {wPmin['phi_mincrit_deg']:.2f}° |

### 3b) Best pair — max (R²_A+R²_B)/2

| | |
|---|---|
| R²_A / R²_B | **{wPmean['R2_A_meancrit']:.4f} / {wPmean['R2_B_meancrit']:.4f}** |
| R²_mean | **{wPmean['R2_mean']:.4f}** |
| Mode | `{wPmean['pair_mean_mode']}` |
| Profile | `{wPmean['profile']}` |
| D/L, f_u | {wPmean['x']:.4f}, {wPmean.get('f_u')} |
| α, d/D | {wPmean['alpha']:.3f}, {wPmean['eps']:.3f} |
| AF | `{wPmean['af_mode']}` |
| κ | **{wPmean['kappa']:.4f}** |
| Amp_A, Amp_B, φ | {wPmean['Amp_A_meancrit']:+.4f}, {wPmean['Amp_B_meancrit']:+.4f}, {wPmean['phi_meancrit_deg']:.2f}° |

## Honest ceiling vs Chris ~0.74

Under **ρ(D)=0.3 ρ_crit** + ST+Gauss H_obs + Amp×sin phase alignment + digitised Fig.6:

- Joint balanced (min criterion) ceiling in this scan is ~**{wPmin['R2_min']:.3f}** (baseline was 0.685).
- Solo A can reach ~**{wA['R2_A_solo']:.3f}**; solo B ~**{wB['R2_B_solo']:.3f}** with independent phases — that is **not** a joint claim.
- Chris sheet ~0.74 likely uses a different AF bookkeeping / κ=0.15 mean-density normalisation / freer dipole; those are **not** the same constraint set.
- Actual κ here remains **~{wPmin['kappa']:.3f}**, not 0.15.

## Method

1. Coarse grid: exponential + two-scale (+ few Plummer), AF ∈ {{hratio, invsq}}, α, d/D.
2. Refine around top seeds (finer x, f_u, α, eps + denser Amp/φ).
3. Pair fits: shared-|Amp| opposite-sign **and** independent Amps with shared φ; take better.
4. Solo fits: independent Amp + φ per series.

## Regenerable

```bash
/workspace/cosmology/.venv/bin/python \\
  /workspace/cosmology/density_local_match/r2_maximize/run_r2_maximize.py
```
"""
    path.write_text(text)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mig = load_migkas()
    print("=== R² maximise under ρ(D)=Ω_m ρ_crit ===", flush=True)

    # Stage 1: coarse
    j1 = stage1_jobs()
    print(f"Stage 1: {len(j1)} physics points (coarse Amp/φ)", flush=True)
    rows1 = run_jobs(j1, mig, refine_phi=False, label="S1")

    # Seeds for refine: top by R2_min, R2_mean, R2_A_solo, R2_B_solo
    seeds = []
    for key in ("R2_min", "R2_mean", "R2_A_solo", "R2_B_solo"):
        ranked = sorted(rows1, key=lambda r: r[key], reverse=True)
        seeds.extend(ranked[:8])
    # unique-ish
    uniq = []
    seen = set()
    for s in seeds:
        k = (s["profile"], round(s["x"], 3), s.get("f_u"), s["af_mode"], round(s["alpha"], 3), round(s["eps"], 3))
        if k not in seen:
            seen.add(k)
            uniq.append(s)
    seeds = uniq[:20]

    j2 = refine_jobs_around(seeds, n_each=8)
    print(f"Stage 2: {len(j2)} refine points (dense Amp/φ)", flush=True)
    rows2 = run_jobs(j2, mig, refine_phi=True, label="S2")

    all_rows = rows1 + rows2

    # Extra tight refine on absolute bests
    best_min_seed = max(all_rows, key=lambda r: (r["R2_min"], r["R2_sum"]))
    best_mean_seed = max(all_rows, key=lambda r: (r["R2_sum"], r["R2_min"]))
    best_a_seed = max(all_rows, key=lambda r: r["R2_A_solo"])
    best_b_seed = max(all_rows, key=lambda r: r["R2_B_solo"])
    j3 = refine_jobs_around([best_min_seed, best_mean_seed, best_a_seed, best_b_seed], n_each=4)
    print(f"Stage 3: {len(j3)} tight refine", flush=True)
    rows3 = run_jobs(j3, mig, refine_phi=True, label="S3")
    all_rows = all_rows + rows3

    # Stage 4: scipy continuous polish around top seeds
    print("Stage 4: scipy polish around top seeds", flush=True)
    try:
        from scipy.optimize import differential_evolution

        def polish_pair_min(seed):
            fam = "twoscale" if seed.get("f_u") is not None or str(seed["profile"]).startswith("twoscale") else (
                "plummer" if str(seed["profile"]).startswith("plummer") else "exponential"
            )
            af = seed["af_mode"]
            if fam == "twoscale":
                bounds = [
                    (max(0.5, seed["x"] - 0.8), seed["x"] + 0.8),
                    (max(0.5, seed["f_u"] - 0.08), min(0.995, seed["f_u"] + 0.08)),
                    (max(0.3, seed["alpha"] - 0.5), seed["alpha"] + 0.5),
                    (max(0.03, seed["eps"] - 0.06), min(0.3, seed["eps"] + 0.06)),
                ]
                def obj(v):
                    x, fu, a, e = map(float, v)
                    try:
                        pmeta = make_pmeta("twoscale", x, f_u=fu)
                        row, _, _, _ = physics_row(pmeta, af, a, e, mig, refine_phi=True)
                        return -row["R2_min"]
                    except Exception:
                        return 0.0
            else:
                bounds = [
                    (max(0.5, seed["x"] - 1.0), seed["x"] + 1.0),
                    (max(0.3, seed["alpha"] - 0.5), seed["alpha"] + 0.5),
                    (max(0.03, seed["eps"] - 0.06), min(0.3, seed["eps"] + 0.06)),
                ]
                def obj(v):
                    x, a, e = map(float, v)
                    try:
                        pmeta = make_pmeta(fam, x)
                        row, _, _, _ = physics_row(pmeta, af, a, e, mig, refine_phi=True)
                        return -row["R2_min"]
                    except Exception:
                        return 0.0
            res = differential_evolution(obj, bounds, maxiter=12, popsize=8, seed=42, polish=True, workers=1)
            v = res.x
            if fam == "twoscale":
                pmeta = make_pmeta("twoscale", float(v[0]), f_u=float(v[1]))
                row, _, _, _ = physics_row(pmeta, af, float(v[2]), float(v[3]), mig, refine_phi=True)
            else:
                pmeta = make_pmeta(fam, float(v[0]))
                row, _, _, _ = physics_row(pmeta, af, float(v[1]), float(v[2]), mig, refine_phi=True)
            return row

        def polish_solo(seed, which="A"):
            fam = "twoscale" if seed.get("f_u") is not None or str(seed["profile"]).startswith("twoscale") else (
                "plummer" if str(seed["profile"]).startswith("plummer") else "exponential"
            )
            af = seed["af_mode"]
            key = "R2_A_solo" if which == "A" else "R2_B_solo"
            if fam == "twoscale":
                bounds = [
                    (max(0.5, seed["x"] - 0.8), seed["x"] + 0.8),
                    (max(0.5, seed["f_u"] - 0.08), min(0.995, seed["f_u"] + 0.08)),
                    (max(0.3, seed["alpha"] - 0.6), seed["alpha"] + 0.6),
                    (max(0.03, seed["eps"] - 0.06), min(0.3, seed["eps"] + 0.06)),
                ]
                def obj(v):
                    x, fu, a, e = map(float, v)
                    try:
                        pmeta = make_pmeta("twoscale", x, f_u=fu)
                        row, _, _, _ = physics_row(pmeta, af, a, e, mig, refine_phi=True)
                        return -row[key]
                    except Exception:
                        return 0.0
            else:
                bounds = [
                    (max(0.5, seed["x"] - 1.2), seed["x"] + 1.2),
                    (max(0.3, seed["alpha"] - 0.6), seed["alpha"] + 0.6),
                    (max(0.03, seed["eps"] - 0.06), min(0.3, seed["eps"] + 0.06)),
                ]
                def obj(v):
                    x, a, e = map(float, v)
                    try:
                        pmeta = make_pmeta(fam, x)
                        row, _, _, _ = physics_row(pmeta, af, a, e, mig, refine_phi=True)
                        return -row[key]
                    except Exception:
                        return 0.0
            res = differential_evolution(obj, bounds, maxiter=12, popsize=8, seed=43, polish=True, workers=1)
            v = res.x
            if fam == "twoscale":
                pmeta = make_pmeta("twoscale", float(v[0]), f_u=float(v[1]))
                row, _, _, _ = physics_row(pmeta, af, float(v[2]), float(v[3]), mig, refine_phi=True)
            else:
                pmeta = make_pmeta(fam, float(v[0]))
                row, _, _, _ = physics_row(pmeta, af, float(v[1]), float(v[2]), mig, refine_phi=True)
            return row

        seed_min = max(all_rows, key=lambda r: (r["R2_min"], r["R2_sum"]))
        seed_mean = max(all_rows, key=lambda r: (r["R2_sum"], r["R2_min"]))
        seed_a = max(all_rows, key=lambda r: r["R2_A_solo"])
        seed_b = max(all_rows, key=lambda r: r["R2_B_solo"])
        for label, fn in [
            ("pair_min", lambda: polish_pair_min(seed_min)),
            ("pair_mean_seed_as_min", lambda: polish_pair_min(seed_mean)),
            ("solo_A", lambda: polish_solo(seed_a, "A")),
            ("solo_B", lambda: polish_solo(seed_b, "B")),
        ]:
            try:
                row = fn()
                all_rows.append(row)
                print(f"  polish {label}: R2_min={row['R2_min']:.4f} RA_solo={row['R2_A_solo']:.4f} RB_solo={row['R2_B_solo']:.4f} κ={row['kappa']:.4f}", flush=True)
            except Exception as ex:
                print(f"  polish {label} failed: {ex}", flush=True)
    except Exception as ex:
        print(f"Stage 4 skipped: {ex}", flush=True)

    best_A = max(all_rows, key=lambda r: r["R2_A_solo"])
    best_B = max(all_rows, key=lambda r: r["R2_B_solo"])
    best_pair_min = max(all_rows, key=lambda r: (r["R2_min"], r["R2_sum"]))
    best_pair_mean = max(all_rows, key=lambda r: (r["R2_sum"], r["R2_min"]))

    top10_joint = sorted(all_rows, key=lambda r: (r["R2_min"], r["R2_sum"]), reverse=True)[:10]
    top5_A = sorted(all_rows, key=lambda r: r["R2_A_solo"], reverse=True)[:5]
    top5_B = sorted(all_rows, key=lambda r: r["R2_B_solo"], reverse=True)[:5]
    top10_mean = sorted(all_rows, key=lambda r: (r["R2_sum"], r["R2_min"]), reverse=True)[:10]

    summary = {
        "constraint": "rho(D)=Omega_m*rho_crit, Omega_m=0.3 (within 1%)",
        "baseline_to_beat": {"R2_A": 0.6852, "R2_B": 0.6855, "kappa": 0.2306},
        "chris_sheet_ballpark": {"R2_A": 0.74, "R2_B": 0.73},
        "n_rows": len(all_rows),
        "best_A_only": slim(best_A),
        "best_B_only": slim(best_B),
        "best_pair_min": slim(best_pair_min),
        "best_pair_mean": slim(best_pair_mean),
        "top10_joint_mincrit": [slim(r) for r in top10_joint],
        "top10_joint_meancrit": [slim(r) for r in top10_mean],
        "top5_A_only": [slim(r) for r in top5_A],
        "top5_B_only": [slim(r) for r in top5_B],
        "script": str(HERE / "run_r2_maximize.py"),
        "method": (
            "ST+Gauss H_obs; AF=hratio|invsq mean-renorm + Δℓ=90°; "
            "AN∝m/M_enc; Amp×sin phase; ρ(D)=0.3ρ_c; report actual κ"
        ),
    }

    # Plots
    # Pair winner (mincrit)
    pmeta_p, phys_p = rebuild_phys(best_pair_min)
    title_p = (
        rf"Best pair (max min $R^2$): {best_pair_min['profile']}, "
        rf"$D/L={best_pair_min['x']:.2f}$, $\alpha={best_pair_min['alpha']:.2f}$, "
        rf"$d/D={best_pair_min['eps']:.2f}$, AF={best_pair_min['af_mode']}, "
        rf"$\kappa={best_pair_min['kappa']:.3f}$  "
        rf"($R^2_A={best_pair_min['R2_A_mincrit']:.3f}$, $R^2_B={best_pair_min['R2_B_mincrit']:.3f}$)"
    )
    plot_winner(
        mig, pmeta_p, phys_p,
        best_pair_min["Amp_A_mincrit"], best_pair_min["Amp_B_mincrit"],
        best_pair_min["phi_mincrit_rad"],
        best_pair_min["R2_A_mincrit"], best_pair_min["R2_B_mincrit"],
        title_p,
        OUT / "migkas_best_pair.png",
        OUT / "migkas_best_pair_scatter.png",
    )

    # A-only plot (show A with solo params; B panel shows B at that same physics with its solo — honest dual panel)
    pmeta_a, phys_a = rebuild_phys(best_A)
    # For A-only PNG: A uses solo A; B uses that physics' own R2_B_solo for context
    plot_solo_combined(
        mig, pmeta_a, phys_a,
        best_A["Amp_A_solo"], best_A["phi_A_solo_rad"],
        best_A["Amp_B_solo"], best_A["phi_B_solo_rad"],
        best_A["R2_A_solo"], best_A["R2_B_solo"],
        (
            f"Best A-only: {best_A['profile']}, D/L={best_A['x']:.2f}, "
            f"α={best_A['alpha']:.2f}, d/D={best_A['eps']:.2f}, "
            f"AF={best_A['af_mode']}, κ={best_A['kappa']:.3f}  "
            f"(R²_A={best_A['R2_A_solo']:.3f} solo; B panel = same phys solo-B)"
        ),
        OUT / "migkas_best_A_only.png",
    )

    pmeta_b, phys_b = rebuild_phys(best_B)
    plot_solo_combined(
        mig, pmeta_b, phys_b,
        best_B["Amp_A_solo"], best_B["phi_A_solo_rad"],
        best_B["Amp_B_solo"], best_B["phi_B_solo_rad"],
        best_B["R2_A_solo"], best_B["R2_B_solo"],
        (
            rf"Best B-only: {best_B['profile']}, $D/L={best_B['x']:.2f}$, "
            rf"$\alpha={best_B['alpha']:.2f}$, $d/D={best_B['eps']:.2f}$, "
            rf"AF={best_B['af_mode']}, $\kappa={best_B['kappa']:.3f}$  "
            rf"($R^2_B={best_B['R2_B_solo']:.3f}$ solo; A panel = same phys solo-A)"
        ),
        OUT / "migkas_best_B_only.png",
    )

    # Also write copies to parent folder as requested names
    for name in ("migkas_best_A_only.png", "migkas_best_B_only.png", "migkas_best_pair.png"):
        src = OUT / name
        dst = PARENT / name
        dst.write_bytes(src.read_bytes())

    summary_path = OUT / "r2_maximize_summary.json"
    # also parent
    summary["beat_baseline_joint"] = bool(best_pair_min["R2_min"] > 0.6852 + 1e-4)
    summary["ceiling_comment"] = (
        f"Joint R2_min={best_pair_min['R2_min']:.4f} vs baseline 0.685; "
        f"solo A={best_A['R2_A_solo']:.4f}, solo B={best_B['R2_B_solo']:.4f}; "
        f"κ_pair={best_pair_min['kappa']:.4f} (not 0.15)."
    )
    summary_path.write_text(json.dumps(summary, indent=2))
    (PARENT / "r2_maximize_summary.json").write_text(json.dumps(summary, indent=2))

    note_path = OUT / "r2_maximize_NOTE.md"
    write_note(summary, note_path)
    (PARENT / "r2_maximize_NOTE.md").write_text(note_path.read_text())

    # Full rows dump (slim)
    (OUT / "r2_maximize_all_rows.json").write_text(
        json.dumps({"n": len(all_rows), "rows": [slim(r) for r in all_rows]}, indent=2)
    )

    print("\n=== WINNERS ===", flush=True)
    print(f"A-only:  R2_A={best_A['R2_A_solo']:.4f}  κ={best_A['kappa']:.4f}  {best_A['profile']} x={best_A['x']:.3f} α={best_A['alpha']} eps={best_A['eps']} af={best_A['af_mode']}", flush=True)
    print(f"B-only:  R2_B={best_B['R2_B_solo']:.4f}  κ={best_B['kappa']:.4f}  {best_B['profile']} x={best_B['x']:.3f} α={best_B['alpha']} eps={best_B['eps']} af={best_B['af_mode']}", flush=True)
    print(f"Pair min: R2_A={best_pair_min['R2_A_mincrit']:.4f} R2_B={best_pair_min['R2_B_mincrit']:.4f} min={best_pair_min['R2_min']:.4f} κ={best_pair_min['kappa']:.4f}", flush=True)
    print(f"Pair mean: R2_A={best_pair_mean['R2_A_meancrit']:.4f} R2_B={best_pair_mean['R2_B_meancrit']:.4f} mean={best_pair_mean['R2_mean']:.4f} κ={best_pair_mean['kappa']:.4f}", flush=True)
    print(f"Wrote under {OUT} and copies in {PARENT}", flush=True)


if __name__ == "__main__":
    main()

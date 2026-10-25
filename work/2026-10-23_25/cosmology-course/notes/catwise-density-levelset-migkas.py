"""CatWISE δ=3% level-set of smooth CoE density laws + Migkas A/A_all refit.

ONLY fixed constraint: δ = [ρ(D−d)−ρ(D+d)]/ρ(D) = 0.031 at primary
d = vision R ≈ 14.3 Gpc.  ρ(D) = Ω_m ρ_crit with Ω_m=0.3.  Smooth entire
way; ρ(0) finite (flag if above stellar-BH density; no floor).

Laws (all δ-locked):
  1. exp   η=1:  ρ=ρ₀ e^{−r/L}                 (L locked; baseline CatWISE-L)
  2. pl:         ρ=ρ₀ (1+r/L)^{−n}              (free n; L from δ)
  3. cored:      ρ=ρ₀/(1+(r/L)²)                (L from δ)
  4. gauss:      ρ=ρ₀ e^{−½(r/L)²}              (L from δ)
  5. soft:       ρ=ρ₀/(1+(r/L)^α)^β             (β=1 fixed; free α; L from δ)
  6. eta (opt.): ρ=ρ₀ e^{−(r/L)^η}  only if interior η≥0.25 improves; skip η≈0.05 edge

Migkas: geometry-only A/A_all = a_f(ℓ−ℓ₀)/⟨a_f⟩, Ob=0, same WLS.
Per law: (a) free ℓ₀  and  (b) ℓ₀ fixed at 15.2° (prefer nearer ptshell).
D freed in Migkas; density-plot baseline D_Earth=58.2 Gpc.

Uses prof_core (kinds exp, pl, gauss, cored, eta, soft).
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
DELTA_USE = 0.031
D_EARTH = 58.2
RV = float(RKM / GPC)          # ≈ 14.3005 Gpc
D_LOCK = RV                    # primary CatWISE lock distance
L0_PTSHELL = 15.2              # preferred ℓ₀ (ptshell ~15.16°)

OMEGA_M = 0.3
H0 = 70.0

CHI2_NULL = 45.13662301218721
N_BINS = 71
PTSHELL_WLS = dict(
    D=58.19699837268826, R2_unweighted=0.5456268309812682,
    R2_weighted=0.5807290671992491, chi2=18.7465256529584,
    l0=15.1629523788697, n_free=5,
)
CATWISE_EXP_D14 = dict(
    D_Gpc=45.99056277928336, L_Gpc=922.6487386323013,
    R2_unweighted=0.716036731924234, chi2=10.572265211558566,
    l0=20.9229393549521, n_free=3,
)
CATWISE_ETA_EDGE = dict(
    D_Gpc=34.96245979193498, eta=0.0500009459948613,
    R2_unweighted=0.745344809732702, chi2=9.278230532227603,
    l0=19.90918750441918, n_free=4,
)
FREE_EXP_WLS = dict(
    D_Gpc=27.95388837552518, L_Gpc=6.072715534614713,
    R2_unweighted=0.7589233227942257, chi2=8.47533016653482, n_free=4,
)

Y = MA
SIG = 0.5 * (MLO + MHI)
W = 1.0 / SIG ** 2
LDEG = np.loadtxt(
    "/workspace/cosmology/migkas/fig6_digitized.csv",
    delimiter=",", skiprows=1, usecols=0,
)

# Stellar-mass BH mean density reference (flag only; no floor).
RHO_BH_MSUN_MPC3 = 1.8e40

# Soft default: β fixed at 1; free α under δ-lock.
SOFT_BETA = 1.0
ETA_FLOOR = 0.25   # skip η-edge; only report if η≥this improves


def rho_crit_Msun_Mpc3(H0_kms: float = H0) -> float:
    return 2.77536627e11 * (H0_kms / 100.0) ** 2


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
# Density laws: shape s(r), δ(D,L,d;shape), L-from-δ, slope, ρ(0)
# ---------------------------------------------------------------------------
def _safe_brentq(f, lo, hi, label=""):
    flo, fhi = f(lo), f(hi)
    if not (np.isfinite(flo) and np.isfinite(fhi)):
        raise RuntimeError(f"non-finite bracket {label}: f({lo})={flo}, f({hi})={fhi}")
    if flo == 0.0:
        return float(lo)
    if fhi == 0.0:
        return float(hi)
    if (flo > 0) == (fhi > 0):
        raise RuntimeError(f"no sign change {label}: f({lo})={flo}, f({hi})={fhi}")
    return float(brentq(f, lo, hi, xtol=1e-12))


def shape_rho(kind, r, L, shape=None):
    """Un-normalised ρ(r); shape dict holds n / alpha / beta / eta."""
    shape = shape or {}
    x = np.asarray(r, dtype=float) / L
    if kind == "exp":
        return np.exp(-x)
    if kind == "eta":
        eta = shape["eta"]
        return np.exp(-(np.abs(x) ** eta))
    if kind == "pl":
        return (1.0 + x) ** (-shape["n"])
    if kind == "cored":
        return 1.0 / (1.0 + x * x)
    if kind == "gauss":
        return np.exp(-0.5 * x * x)
    if kind == "soft":
        a, b = shape["alpha"], shape["beta"]
        return (1.0 + np.abs(x) ** a) ** (-b)
    raise ValueError(kind)


def delta_of(kind, D, L, d, shape=None):
    if L <= 0 or D <= d or d <= 0:
        return np.nan
    sD = float(shape_rho(kind, D, L, shape))
    if sD <= 0 or not np.isfinite(sD):
        return np.nan
    return float(
        (shape_rho(kind, D - d, L, shape) - shape_rho(kind, D + d, L, shape)) / sD
    )


def L_from_delta(kind, D, d, shape=None, delta: float = DELTA_USE) -> float:
    """Solve δ(D,L,d;shape)=delta for L>0."""
    shape = shape or {}
    if kind == "exp" or (kind == "eta" and abs(shape.get("eta", 1.0) - 1.0) < 1e-12):
        return float(d / np.arcsinh(delta / 2.0))

    def f(logL):
        L = np.exp(logL)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                return delta_of(kind, D, L, d, shape) - delta

    # Bracket on log L: large L → δ→0; small L → large δ (usually)
    grid = np.linspace(np.log(1e-4), np.log(1e6), 2500)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        vals = np.array([f(g) for g in grid])
    ok = np.isfinite(vals)
    for i in range(len(grid) - 1):
        if not (ok[i] and ok[i + 1]):
            continue
        a, b = vals[i], vals[i + 1]
        if a == 0.0 or b == 0.0 or (a > 0) != (b > 0):
            logL = _safe_brentq(f, float(grid[i]), float(grid[i + 1]),
                                label=f"{kind} D={D}")
            L = float(np.exp(logL))
            if L > 0 and np.isfinite(L):
                return L
    raise RuntimeError(f"no L root for kind={kind} D={D} d={d} shape={shape}")


def slope_ln_at_D(kind, D, L, shape=None) -> float:
    """(ρ'/ρ)_D by central difference (analytic where cheap)."""
    shape = shape or {}
    if kind == "exp":
        return -1.0 / L
    if kind == "eta":
        eta = shape["eta"]
        return float(-(eta / L) * (D / L) ** (eta - 1.0))
    if kind == "pl":
        n = shape["n"]
        return float(-n / (L + D))
    if kind == "cored":
        return float(-2.0 * D / (L * L + D * D))
    if kind == "gauss":
        return float(-D / (L * L))
    if kind == "soft":
        a, b = shape["alpha"], shape["beta"]
        # d/dr log(1+(r/L)^a)^{-b} = -b * a (r/L)^{a-1} (1/L) / (1+(r/L)^a)
        x = D / L
        return float(-b * a * (x ** (a - 1.0)) / L / (1.0 + x ** a))
    # numerical fallback
    eps = 1e-4 * D
    s1 = float(shape_rho(kind, D - eps, L, shape))
    s2 = float(shape_rho(kind, D + eps, L, shape))
    s0 = float(shape_rho(kind, D, L, shape))
    return float((s2 - s1) / (2 * eps) / s0)


def rho0_over_crit(kind, D, L, shape=None) -> float:
    """ρ(0)/ρ_crit given ρ(D)=Ω_m ρ_crit."""
    s0 = float(shape_rho(kind, 0.0, L, shape))
    sD = float(shape_rho(kind, D, L, shape))
    return float(OMEGA_M * s0 / sD)


def rho_over_crit(kind, r, D, L, shape=None):
    s = shape_rho(kind, r, L, shape)
    sD = float(shape_rho(kind, D, L, shape))
    return OMEGA_M * s / sD


def density_card(kind, D, d, shape=None, delta: float = DELTA_USE):
    shape = dict(shape or {})
    L = L_from_delta(kind, D, d, shape, delta)
    dcheck = delta_of(kind, D, L, d, shape)
    slope = slope_ln_at_D(kind, D, L, shape)
    r0 = rho0_over_crit(kind, D, L, shape)
    rho_c = rho_crit_Msun_Mpc3()
    below_bh = bool(r0 * rho_c < RHO_BH_MSUN_MPC3)
    # gentleness: |Δ slope| across D±d/2 relative to |slope| — small ⇒ gently varying
    half = 0.5 * d
    s_lo = slope_ln_at_D(kind, max(D - half, 0.01 * D), L, shape)
    s_hi = slope_ln_at_D(kind, D + half, L, shape)
    cliff = float(abs(s_hi - s_lo) / max(abs(slope), 1e-30))
    return dict(
        kind=kind,
        d_Gpc=float(d),
        D_Gpc=float(D),
        L_Gpc=float(L),
        L_over_D=float(L / D),
        shape=shape,
        delta_target=float(delta),
        delta_check=float(dcheck),
        delta_rel_err=float(abs(dcheck - delta) / delta),
        Omega_m=float(OMEGA_M),
        rho_D_over_rho_crit=float(OMEGA_M),
        rho_0_over_rho_crit=float(r0),
        rho_prime_over_rho_at_D=float(slope),
        drho_dr_at_D_over_rho_crit_per_Gpc=float(OMEGA_M * slope),
        cliff_index=cliff,
        rho_0_below_BH_density=below_bh,
        rho_BH_ref_Msun_Mpc3=float(RHO_BH_MSUN_MPC3),
        law=_law_latex(kind, shape),
    )


def _g(v):
    return '?' if v is None else f"{v:.3g}"


def _law_latex(kind, shape):
    shape = shape or {}
    if kind == "exp":
        return r"ρ=ρ₀ e^{−r/L}"
    if kind == "eta":
        return rf"ρ=ρ₀ e^{{−(r/L)^{_g(shape.get('eta'))}}}"
    if kind == "pl":
        return rf"ρ=ρ₀ (1+r/L)^{{−{_g(shape.get('n'))}}}"
    if kind == "cored":
        return r"ρ=ρ₀/(1+(r/L)²)"
    if kind == "gauss":
        return r"ρ=ρ₀ e^{−½(r/L)²}"
    if kind == "soft":
        return (
            rf"ρ=ρ₀/(1+(r/L)^{_g(shape.get('alpha'))})"
            rf"^{{{_g(shape.get('beta'))}}}"
        )
    return kind


# ---------------------------------------------------------------------------
# Migkas sky model via prof_core
# ---------------------------------------------------------------------------
def _pc_args(kind, DR, LD, shape):
    a = dict(DR=float(DR), LD=float(LD))
    shape = shape or {}
    if kind == "pl":
        a["n"] = float(shape["n"])
    elif kind == "eta":
        a["eta"] = float(shape["eta"])
    elif kind == "soft":
        a["alpha"] = float(shape["alpha"])
        a["beta"] = float(shape["beta"])
    return a


def sky_locked(kind, D, mc, d_lock, shape=None):
    """Mean-normalised a_f; L from CatWISE δ lock.  Returns (af, rr, LD, L) or None."""
    shape = dict(shape or {})
    if not (D > max(d_lock, RV) * 1.05 and 1e-4 < mc < 80.0):
        return None
    # shape bounds
    if kind == "pl" and not (1.2 < shape.get("n", 0) < 10.0):
        return None
    if kind == "soft" and not (0.3 < shape.get("alpha", 0) < 8.0):
        return None
    if kind == "eta" and not (ETA_FLOOR <= shape.get("eta", 0) < 8.0):
        return None
    try:
        L = L_from_delta(kind, D, d_lock, shape, DELTA_USE)
    except Exception:
        return None
    if not (np.isfinite(L) and L > 1e-3):
        return None
    LD = L / D
    if LD < 1e-4:
        return None
    DR = D / RV
    a = _pc_args(kind, DR, LD, shape)
    try:
        rr = pc.run(kind if kind != "exp" else "exp", a, float(mc), 0.0)
    except Exception:
        return None
    af = rr["af"]
    if (not np.all(np.isfinite(af))) or rr["vrmin"] <= 0:
        return None
    return af / af.mean(), rr, LD, L


def model_bins(kind, D, mc, l0, d_lock, shape=None):
    out = sky_locked(kind, D, mc, d_lock, shape)
    if out is None:
        return None
    af, rr, LD, L = out
    return shifted(af, l0), rr, LD, L


def chi2_of(kind, D, mc, l0, d_lock, shape=None):
    out = model_bins(kind, D, mc, l0, d_lock, shape)
    if out is None:
        return 800.0
    pred, _, _, _ = out
    return float(np.sum(W * (pred[IMIG] - Y) ** 2))


# ---------------------------------------------------------------------------
# Fit one law: free ℓ₀ and/or fixed ℓ₀=15.2°
# ---------------------------------------------------------------------------
def _pack_fit(kind, D, mc, l0, L, LD, rr, full, d_lock, shape, n_free, mode):
    m = full[IMIG]
    c2 = float(np.sum(W * (m - Y) ** 2))
    dcheck = delta_of(kind, D, L, d_lock, shape)
    r0 = rho0_over_crit(kind, D, L, shape)
    rho_c = rho_crit_Msun_Mpc3()
    below_bh = bool(r0 * rho_c < RHO_BH_MSUN_MPC3)
    dchi2 = CHI2_NULL - c2
    sig, pval = sigma_from_delta(dchi2, n_free)
    half = 0.5 * d_lock
    s0 = slope_ln_at_D(kind, D, L, shape)
    s_lo = slope_ln_at_D(kind, max(D - half, 0.01 * D), L, shape)
    s_hi = slope_ln_at_D(kind, D + half, L, shape)
    cliff = float(abs(s_hi - s_lo) / max(abs(s0), 1e-30))
    out = dict(
        kind=kind,
        mode=mode,
        D_Gpc=float(D),
        DR=float(D / RV),
        L_Gpc_locked=float(L),
        LD=float(LD),
        mc=float(mc),
        Ob=0.0,
        l0=float(l0 % 360.0),
        shape=dict(shape or {}),
        chi2=c2,
        chi2_red=float(c2 / max(Y.size - n_free, 1)),
        R2_unweighted=r2_u(Y, m),
        R2_weighted=r2_w(Y, m, W),
        A_min=float(full.min()),
        A_max=float(full.max()),
        A_amp=float(full.max() - full.min()),
        vrmin=float(rr["vrmin"]),
        Ht=float(rr["Ht"]),
        age_Gyr=float(rr["age"]),
        contrast_Rvis_model=float(rr["contrast"]),
        contrast_analytic_at_d=float(dcheck),
        d_used_for_lock_Gpc=float(d_lock),
        rho_0_over_rho_crit=float(r0),
        rho_0_below_BH_density=below_bh,
        rho_prime_over_rho_at_D=float(s0),
        cliff_index=cliff,
        n_free=int(n_free),
        Delta_chi2_vs_null=float(dchi2),
        Delta_nu_vs_null=int(n_free),
        p_value_nested=float(pval) if pval is not None else None,
        sigma_nested_one_sided=float(sig),
        sigma_sqrt_Delta_chi2=float(np.sqrt(max(dchi2, 0.0))),
        chi2_null=float(CHI2_NULL),
        N_bins=int(N_BINS),
        law=_law_latex(kind, shape),
        fixed_L_rule=(
            "L solved from δ=[ρ(D−d)−ρ(D+d)]/ρ(D)=0.031 for each trial "
            f"(D, shape); d={d_lock:.4g} Gpc"
        ),
    )
    return out


def fit_law(kind, d_lock: float, l0_mode: str = "free",
            shape_seed: dict | None = None):
    """l0_mode: 'free' | 'fixed15.2'.  Shape params free where applicable."""
    free_shape = kind in ("pl", "soft", "eta")
    fix_l0 = (l0_mode == "fixed15.2")
    l0_fixed = L0_PTSHELL

    # Parameter vector layout:
    #   always: log10(D), log10(mc)
    #   + l0 if free
    #   + shape param if free_shape (log10(n) | log10(alpha) | log10(eta))
    # Compact start grid (paired D/mc/l0/shape) to keep runtime manageable.
    starts_Dmc = [
        (58.2, 0.2), (46.0, 4.5), (40.0, 0.5), (30.0, 1.0),
        (50.0, 2.0), (35.0, 0.8), (55.0, 1.5), (28.0, 0.2),
        (70.0, 0.1), (42.0, 5.0),
    ]
    l0_seeds = [15.2, 14.0, 16.0, 18.0, 20.0, 12.0, 21.0, 13.0, 15.0, 19.0]
    shape_seeds = {
        "pl": [2.0, 3.0, 4.0, 1.5, 6.0],
        "soft": [1.0, 1.5, 2.0, 0.8, 2.5],
        "eta": [0.5, 0.7, 1.0, 0.35, 1.5],
    }

    starts = []
    if free_shape:
        # ~10 paired starts covering D/mc/l0/shape (not full Cartesian)
        sh_list = shape_seeds[kind]
        for i, (D0, mc0) in enumerate(starts_Dmc):
            l0_0 = l0_fixed if fix_l0 else l0_seeds[i % len(l0_seeds)]
            s0 = sh_list[i % len(sh_list)]
            u = [np.log10(D0), np.log10(mc0)]
            if not fix_l0:
                u.append(l0_0)
            u.append(np.log10(s0))
            starts.append(np.array(u, dtype=float))
        # a few extra shape probes at the CatWISE-exp-like (D,mc)
        for s0 in sh_list:
            u = [np.log10(46.0), np.log10(4.5)]
            if not fix_l0:
                u.append(15.2 if fix_l0 else 20.0)
            u.append(np.log10(s0))
            starts.append(np.array(u, dtype=float))
    else:
        for i, (D0, mc0) in enumerate(starts_Dmc):
            l0_0 = l0_fixed if fix_l0 else l0_seeds[i % len(l0_seeds)]
            u = [np.log10(D0), np.log10(mc0)]
            if not fix_l0:
                u.append(l0_0)
            starts.append(np.array(u, dtype=float))

    # seed from known CatWISE-exp / ptshell-ish
    if kind == "exp" and not free_shape:
        extra = [np.array([np.log10(45.99), np.log10(4.53)]
                          + ([] if fix_l0 else [20.92]))]
        starts = extra + starts

    def unpack(u):
        D = 10 ** u[0]
        mc = 10 ** u[1]
        idx = 2
        if fix_l0:
            l0 = l0_fixed
        else:
            l0 = u[idx]
            idx += 1
        shape = dict(shape_seed or {})
        if kind == "pl":
            shape = dict(n=float(10 ** u[idx]))
        elif kind == "soft":
            shape = dict(alpha=float(10 ** u[idx]), beta=SOFT_BETA)
        elif kind == "eta":
            shape = dict(eta=float(10 ** u[idx]))
        return D, mc, l0, shape

    def objective(u):
        D, mc, l0, shape = unpack(u)
        return chi2_of(kind, D, mc, l0, d_lock, shape)

    n_free = 2 + (0 if fix_l0 else 1) + (1 if free_shape else 0)
    tag = f"{kind}/{l0_mode}"
    best_u, best_c = None, np.inf
    for i, u0 in enumerate(starts):
        D0, mc0, l00, sh0 = unpack(u0)
        shs = ",".join(f"{k}={v:.3g}" for k, v in sh0.items()) or "-"
        print(
            f"  [{tag}] start {i+1}/{len(starts)}: "
            f"D={D0:.1f} mc={mc0:.4g} l0={l00:.1f} shape={shs}",
            flush=True,
        )
        res = minimize(
            objective, u0, method="Nelder-Mead",
            options=dict(maxiter=260, xatol=1e-4, fatol=1e-4, adaptive=True),
        )
        print(f"    -> chi2={res.fun:.4f} nit={res.nit}", flush=True)
        if res.fun < best_c:
            best_c, best_u = float(res.fun), res.x.copy()

    if best_u is None:
        raise RuntimeError(f"no viable fit for {tag}")
    D, mc, l0, shape = unpack(best_u)
    out = model_bins(kind, D, mc, l0, d_lock, shape)
    if out is None:
        raise RuntimeError(f"best params illegal for {tag}")
    full, rr, LD, L = out
    fitp = _pack_fit(kind, D, mc, l0, L, LD, rr, full, d_lock, shape, n_free,
                     l0_mode)
    fitp["free_params"] = (
        ["D", "mc"] + ([] if fix_l0 else ["l0"])
        + (["n"] if kind == "pl" else [])
        + (["alpha"] if kind == "soft" else [])
        + (["eta"] if kind == "eta" else [])
    )
    return fitp, full


# ---------------------------------------------------------------------------
# Density overlay plot (level set at D=58.2)
# ---------------------------------------------------------------------------
def plot_density_levelset(cards: dict, best_shapes: dict | None = None):
    """Overlay all δ=3% level-set curves; annotate Earth / D±d / slopes."""
    D = D_EARTH
    d = D_LOCK
    r_max = D + 1.55 * d
    r = np.linspace(0.0, r_max, 3000)

    # Display order and styles
    order = ["exp", "pl", "cored", "gauss", "soft"]
    if best_shapes and "eta" in best_shapes:
        order.append("eta")
    colors = {
        "exp": "C0", "pl": "C1", "cored": "C2", "gauss": "C4",
        "soft": "C5", "eta": "C3",
    }
    lstyles = {
        "exp": "-", "pl": "--", "cored": "-.", "gauss": ":",
        "soft": (0, (5, 1, 1, 1)), "eta": (0, (3, 1, 1, 1, 1, 1)),
    }

    fig, ax = plt.subplots(figsize=(14.0, 7.2))
    csv_cols = [r]
    csv_hdr = ["r_Gpc"]
    ymax = 0.0
    params_out = {}

    # Use best-fit shapes for free-shape laws when available; else defaults
    default_shapes = {
        "exp": {},
        "pl": {"n": 3.0},
        "cored": {},
        "gauss": {},
        "soft": {"alpha": 1.5, "beta": SOFT_BETA},
        "eta": {"eta": 0.5},
    }
    shapes_use = dict(default_shapes)
    if best_shapes:
        for k, sh in best_shapes.items():
            shapes_use[k] = sh

    primary_kind = None
    for kind in order:
        shape = shapes_use[kind]
        try:
            card = density_card(kind, D, d, shape)
        except Exception as e:
            print(f"  skip density {kind}: {e}", flush=True)
            continue
        params_out[kind] = card
        L = card["L_Gpc"]
        y = rho_over_crit(kind, r, D, L, shape)
        ymax = max(ymax, float(np.nanmax(y[np.isfinite(y)])))
        shs = ",".join(f"{kk}={vv:.3g}" for kk, vv in shape.items())
        gentle = "gentle" if card["cliff_index"] < 0.5 else (
            "mild" if card["cliff_index"] < 1.5 else "steep/cliffy"
        )
        lab = (
            rf"{kind}: $L={L:.1f}$ Gpc"
            + (rf" ({shs})" if shs else "")
            + rf", $\delta={card['delta_check']:.4f}$, cliff={card['cliff_index']:.2f} ({gentle})"
        )
        ax.plot(r, y, color=colors[kind], ls=lstyles[kind], lw=2.4, label=lab)
        csv_cols.append(y)
        csv_hdr.append(
            f"rho_over_rhocrit_{kind}_L{L:.4g}"
            + (f"_{shs.replace(',', '_').replace('=', '')}" if shs else "")
        )
        if primary_kind is None:
            primary_kind = kind

    # Annotate Earth / D±d using exp card (same δ by construction)
    card0 = params_out.get("exp") or list(params_out.values())[0]
    L0 = card0["L_Gpc"]
    kind0 = card0["kind"]
    shape0 = card0["shape"]
    slope = card0["drho_dr_at_D_over_rho_crit_per_Gpc"]

    ax.axvline(D, color="k", lw=1.8, ls="--",
               label=rf"Earth at $D={D:.1f}$ Gpc (CoE→Earth)")
    ax.axvline(D - d, color="0.35", lw=1.1, ls=":")
    ax.axvline(D + d, color="0.35", lw=1.1, ls=":",
               label=rf"$D\pm d$ with $d={d:.2f}$ Gpc (vision $R$)")
    ax.axhline(OMEGA_M, color="0.45", lw=0.9, alpha=0.8)
    ax.plot([D], [OMEGA_M], "o", color="k", ms=8, zorder=5)

    # Tangent of exp (representative local slope under δ lock — nearly universal)
    half = 1.2 * d
    rt = np.linspace(max(0.0, D - half), D + half, 200)
    yt = OMEGA_M + slope * (rt - D)
    ax.plot(
        rt, yt, color="0.25", lw=1.6, ls="-.",
        label=(
            rf"Tangent (exp): $(d\rho/dr)_D={slope:.3e}\,\rho_{{\rm crit}}/$Gpc"
            rf"; $\rho(D)/\rho_{{\rm crit}}=\Omega_m={OMEGA_M:g}$"
        ),
    )

    rho_m = float(rho_over_crit(kind0, D - d, D, L0, shape0))
    rho_p = float(rho_over_crit(kind0, D + d, D, L0, shape0))
    ax.plot([D - d], [rho_m], "s", color="0.2", ms=7, zorder=5)
    ax.plot([D + d], [rho_p], "s", color="0.2", ms=7, zorder=5)

    # Bottom-band notes (not on curves)
    ax.annotate(
        rf"$(d\rho/dr)_{{D}}\approx{slope:.3e}\,\rho_{{\rm crit}}/$Gpc"
        f"\n"
        rf"(local slope ≈ universal under $\delta$ lock)",
        xy=(D, OMEGA_M),
        xytext=(30.0, 0.247),
        fontsize=11, ha="left", va="bottom",
        arrowprops=dict(arrowstyle="->", color="0.2", lw=1.0,
                        connectionstyle="arc3,rad=0.3"),
        bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.92),
        zorder=6,
    )
    ax.annotate(
        rf"$\rho(D-d)/\rho_{{\rm crit}}\approx{rho_m:.4f}$",
        xy=(D - d, rho_m),
        xytext=(D - d + 2.0, 0.272),
        fontsize=10, ha="center",
        arrowprops=dict(arrowstyle="->", color="0.3", lw=0.9),
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                  edgecolor="0.4", alpha=0.92),
        zorder=6,
    )
    ax.annotate(
        rf"$\rho(D+d)/\rho_{{\rm crit}}\approx{rho_p:.4f}$",
        xy=(D + d, rho_p),
        xytext=(min(D + d - 6.0, r_max - 14.0), 0.262),
        fontsize=10, ha="left", va="bottom",
        arrowprops=dict(arrowstyle="->", color="0.3", lw=0.9),
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                  edgecolor="0.4", alpha=0.92),
        zorder=6,
    )

    # ρ(0) callouts for each (compact, bottom-left)
    r0_lines = []
    for kind in order:
        if kind not in params_out:
            continue
        c = params_out[kind]
        flag = "OK" if c["rho_0_below_BH_density"] else "FLAG>BH"
        r0_lines.append(f"{kind}: ρ(0)/ρc={c['rho_0_over_rho_crit']:.3f} ({flag})")
    ax.text(
        0.01, 0.015, "ρ(0) finite:\n" + "\n".join(r0_lines),
        transform=ax.transAxes, fontsize=9, va="bottom", ha="left",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="lightyellow", alpha=0.92),
        family="monospace",
    )

    ax.set_xlabel(
        r"Distance from CoE (Centre of Expansion), $r$ (Gpc)", fontsize=13,
    )
    ax.set_ylabel(r"Mass density $\rho(r)/\rho_{\mathrm{crit}}$", fontsize=13)
    ax.set_title(
        r"CatWISE $\delta=0.031$ level set of smooth CoE densities "
        r"(no cliffs / shell edges)"
        "\n"
        rf"$d=R_{{\rm vis}}={d:.2f}$ Gpc locks $L$ (and free-shape params "
        rf"shown at best Migkas / default); "
        rf"$\rho(D)=\Omega_m\rho_{{\rm crit}}$, $\Omega_m={OMEGA_M:g}$, "
        rf"$D={D:.1f}$ Gpc",
        fontsize=13,
    )
    ax.set_xlim(0, r_max)
    # Reference: rejected η≈0.05 edge (cliff near CoE then flat) — faint, for contrast
    try:
        sh_e = {"eta": 0.05}
        Le = L_from_delta("eta", D, d, sh_e)
        ye = rho_over_crit("eta", r, D, Le, sh_e)
        de = delta_of("eta", D, Le, d, sh_e)
        ax.plot(r, ye, color="0.6", lw=1.6, ls="-", alpha=0.8, zorder=1,
                label=(rf"REJECTED η≈0.05 edge (cliff-then-flat): $\delta={de:.4f}$, "
                       rf"$\rho(0)/\rho_{{\rm crit}}={float(ye[0]):.2f}$ (off-scale)"))
    except Exception as e:
        print("  edge ref skipped:", e)
    ax.set_ylim(0.24, 0.42)
    ax.tick_params(labelsize=12)
    ax.grid(alpha=0.25)

    fig.subplots_adjust(left=0.07, right=0.52, top=0.88, bottom=0.11)
    ax.legend(
        fontsize=14, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        ncol=1, framealpha=0.95, borderaxespad=0.0,
        fancybox=True, edgecolor="0.6",
    )

    png = OUT / "catwise-density-levelset.png"
    fig.savefig(png, dpi=160, bbox_inches="tight")
    plt.close(fig)

    csv = OUT / "catwise-density-levelset.csv"
    np.savetxt(
        csv, np.column_stack(csv_cols), delimiter=",",
        header=",".join(csv_hdr), comments="",
    )
    print("wrote", png)
    print("wrote", csv)
    return params_out


# ---------------------------------------------------------------------------
# Migkas best-curve PNGs
# ---------------------------------------------------------------------------
def plot_migkas(fitp, full, tag: str):
    fig, ax = plt.subplots(figsize=(19.0, 6.0))
    ax.errorbar(
        LDEG, Y, yerr=np.vstack([MLO, MHI]), fmt="o", ms=4.0, color="k",
        ecolor="0.45", elinewidth=0.7, capsize=0, zorder=3,
        label=r"Migkas Fig. 6  $A/A_{\mathrm{all}}$",
    )
    sh = fitp.get("shape") or {}
    shs = ", ".join(rf"${k}={v:.3g}$" for k, v in sh.items())
    ax.plot(
        TH1, full, color="C0", lw=2.0,
        label=(
            rf"{fitp['kind']} ($\delta$ lock): $L={fitp['L_Gpc_locked']:.1f}$ Gpc"
            + (f", {shs}" if shs else "")
            + f"\n"
            + rf"$D={fitp['D_Gpc']:.1f}$ Gpc, $L/D={fitp['LD']:.3f}$, "
            + rf"$m_c={fitp['mc']:.4g}$, $\ell_0={fitp['l0']:.1f}^\circ$"
            + f" [{fitp['mode']}]"
            + f"\n"
            + rf"$R^2={fitp['R2_unweighted']:.3f}$ (unw.), "
            + rf"$\chi^2={fitp['chi2']:.2f}$; "
            + rf"$\sigma_{{\rm null}}={fitp['sigma_nested_one_sided']:.2f}$; "
            + rf"$A_{{\rm amp}}={fitp['A_amp']:.4f}$"
            + f"\n"
            + rf"CatWISE-exp $R^2\approx{CATWISE_EXP_D14['R2_unweighted']:.2f}$; "
            + rf"ptshell $R^2\approx{PTSHELL_WLS['R2_unweighted']:.2f}$"
        ),
    )
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)", fontsize=13)
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$", fontsize=13)
    ax.set_title(
        rf"Migkas geometry-only, CatWISE $\delta={DELTA_USE:g}$ lock — "
        rf"{fitp['kind']} / {fitp['mode']}"
        f"\n"
        rf"Ob=0; free {{{', '.join(fitp['free_params'])}}}; $L=L(D,\mathrm{{shape}},d)$",
        fontsize=14,
    )
    ax.legend(fontsize=14, loc="upper left", bbox_to_anchor=(1.02, 1.0), framealpha=0.95)
    ax.set_xticks(np.arange(0, 361, 45))
    ax.tick_params(labelsize=12)
    ax.grid(alpha=0.25)
    fig.subplots_adjust(left=0.05, right=0.56, top=0.86, bottom=0.12)
    png = OUT / f"catwise-density-levelset-migkas-{tag}.png"
    fig.savefig(png, dpi=160, bbox_inches="tight")
    plt.close(fig)
    csv = OUT / f"catwise-density-levelset-migkas-{tag}.csv"
    np.savetxt(
        csv, np.column_stack([TH1, full]), delimiter=",",
        header="ell_deg,A_over_Aall", comments="",
    )
    print("wrote", png)
    print("wrote", csv)
    return png


def plot_migkas_overlay(fits_curves: dict, mode: str):
    """Overlay best curves for all laws under one ℓ₀ mode."""
    fig, ax = plt.subplots(figsize=(11.0, 5.8))
    ax.errorbar(
        LDEG, Y, yerr=np.vstack([MLO, MHI]), fmt="o", ms=3.5, color="k",
        ecolor="0.45", elinewidth=0.6, capsize=0, zorder=3,
        label=r"Migkas Fig. 6",
    )
    colors = {
        "exp": "C0", "pl": "C1", "cored": "C2", "gauss": "C4",
        "soft": "C5", "eta": "C3",
    }
    for kind, (fitp, full) in fits_curves.items():
        ax.plot(
            TH1, full, color=colors.get(kind, "C0"), lw=1.8,
            label=(
                rf"{kind}: $R^2={fitp['R2_unweighted']:.3f}$, "
                rf"$\chi^2={fitp['chi2']:.2f}$, "
                rf"$\ell_0={fitp['l0']:.1f}^\circ$"
            ),
        )
    ax.set_xlim(0, 360)
    ax.set_xlabel(r"Galactic longitude $\ell$ (deg)", fontsize=13)
    ax.set_ylabel(r"$A/A_{\mathrm{all}}$", fontsize=13)
    ax.set_title(
        rf"Migkas level-set overlay — $\ell_0$ mode: {mode}; "
        rf"CatWISE $\delta={DELTA_USE:g}$ lock"
        f"\n"
        r"geometry-only, Ob=0; each law $L$ locked by $\delta$",
        fontsize=14,
    )
    ax.legend(
        fontsize=14, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        framealpha=0.95,
    )
    ax.set_xticks(np.arange(0, 361, 45))
    ax.tick_params(labelsize=12)
    ax.grid(alpha=0.25)
    fig.subplots_adjust(right=0.72)
    fig.tight_layout()
    png = OUT / f"catwise-density-levelset-migkas-overlay-{mode}.png"
    fig.savefig(png, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", png)
    return png


# ---------------------------------------------------------------------------
# Markdown note
# ---------------------------------------------------------------------------
def write_md(dens_params: dict, fits: dict):
    lines = []
    lines.append("# CatWISE δ=3% density level set + Migkas refit")
    lines.append("")
    lines.append("**Date:** 2026-10-06 (America/Chicago)")
    lines.append("")
    lines.append(
        "**Acronyms:** CatWISE = CatWISE quasar dipole (Secrest et al. 2021); "
        "CoE = centre of expansion; Migkas = Migkas et al. cluster "
        r"$A/A_{\mathrm{all}}$ anisotropy."
    )
    lines.append("")
    lines.append("## Constraint (only fixed)")
    lines.append("")
    lines.append(
        r"$$\delta=\frac{\rho(D-d)-\rho(D+d)}{\rho(D)}=0.031"
        r"\quad\text{at primary }d=R_{\rm vis}\approx 14.3\,\mathrm{Gpc}.$$"
    )
    lines.append("")
    lines.append(
        rf"- $\rho(D)=\Omega_m\rho_{{\rm crit}}$ with $\Omega_m={OMEGA_M}$."
    )
    lines.append(
        f"- Density-plot baseline $D={D_EARTH}$ Gpc; Migkas fit frees $D$."
    )
    lines.append(
        r"- Smooth entire way; $\rho(0)$ finite (flag if above stellar-BH "
        r"density; no floor)."
    )
    lines.append(
        r"- Prefer $\ell_0$ nearer ptshell $\sim 15.2^\circ$ "
        r"(Chris: $19.9^\circ$ too much): report free $\ell_0$ and "
        rf"$\ell_0={L0_PTSHELL}^\circ$ fixed."
    )
    lines.append("")
    lines.append("## Density family (level set at $D=58.2$ Gpc)")
    lines.append("")
    lines.append(
        "| law | formula | shape | $L$ (Gpc) | "
        r"$\delta$ check | $\rho(0)/\rho_c$ | "
        r"$(\rho'/\rho)_D$ | cliff index | gentle? |"
    )
    lines.append("|---|---|---|---:|---:|---:|---:|---:|---|")
    for kind, c in dens_params.items():
        sh = ", ".join(f"{k}={v:.3g}" for k, v in (c.get("shape") or {}).items()) or "—"
        gentle = "yes" if c["cliff_index"] < 0.5 else (
            "mild" if c["cliff_index"] < 1.5 else "cliffy"
        )
        lines.append(
            f"| {kind} | `{c['law']}` | {sh} | {c['L_Gpc']:.4g} | "
            f"{c['delta_check']:.5f} | {c['rho_0_over_rho_crit']:.4g} | "
            f"{c['rho_prime_over_rho_at_D']:.6g} | {c['cliff_index']:.3g} | "
            f"{gentle} |"
        )
    lines.append("")
    lines.append(
        "**Cliff index** $=|(\\rho'/\\rho)_{D+d/2}-(\\rho'/\\rho)_{D-d/2}|$"
        r"$/|(\rho'/\rho)_D|$: small ⇒ slope varies gently through Earth "
        "(no precipice). The η-edge solution (η≈0.05) is cliffy by design "
        "and is **not** treated as winner here."
    )
    lines.append("")

    # Ranking tables
    for mode, title in [
        ("fixed15.2", rf"Ranking under $\ell_0={L0_PTSHELL}^\circ$ fixed"),
        ("free", r"Ranking under free $\ell_0$"),
    ]:
        rows = []
        for key, fp in fits.items():
            if fp.get("mode") != mode:
                continue
            rows.append(fp)
        rows.sort(key=lambda x: x["chi2"])
        lines.append(f"## {title}")
        lines.append("")
        lines.append(
            r"| rank | law | $R^2$ | $\chi^2$ | $n_{\rm free}$ | $\Delta\chi^2$ | "
            r"$\sigma$ null (nested) | $\sqrt{\Delta\chi^2}$ | "
            r"$D$ | $L$ | shape | $m_c$ | $\ell_0$ | cliff | beats 0.72? |"
        )
        lines.append("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---|")
        for i, fp in enumerate(rows, 1):
            sh = ", ".join(
                f"{k}={v:.3g}" for k, v in (fp.get("shape") or {}).items()
            ) or "—"
            beats = "YES" if fp["R2_unweighted"] > 0.72 and fp["l0"] <= 16.0 else (
                "R²>0.72 but ℓ₀>16°" if fp["R2_unweighted"] > 0.72 else "no"
            )
            lines.append(
                f"| {i} | {fp['kind']} | {fp['R2_unweighted']:.4f} | "
                f"{fp['chi2']:.3f} | {fp['n_free']} | {fp['Delta_chi2_vs_null']:.3f} | "
                f"{fp['sigma_nested_one_sided']:.2f} | {fp['sigma_sqrt_Delta_chi2']:.2f} | "
                f"{fp['D_Gpc']:.2f} | {fp['L_Gpc_locked']:.4g} | {sh} | "
                f"{fp['mc']:.4g} | {fp['l0']:.2f} | {fp['cliff_index']:.3g} | "
                f"{beats} |"
            )
        lines.append("")
        if rows:
            best = rows[0]
            lines.append(
                f"**Best under {mode}:** `{best['kind']}` with "
                f"$R^2={best['R2_unweighted']:.4f}$, "
                f"$\\chi^2={best['chi2']:.3f}$, "
                f"$\\ell_0={best['l0']:.2f}^\\circ$, "
                f"cliff={best['cliff_index']:.3g}."
            )
            lines.append("")

    lines.append("## σ vs null for every fit")
    lines.append("")
    lines.append(
        rf"Null: flat $A/A_{{\rm all}}=1$, $\chi^2_{{\rm null}}={CHI2_NULL:.3f}$, $N={N_BINS}$. "
        r"Nested: $\Delta\chi^2=\chi^2_{\rm null}-\chi^2$, $p=P(\chi^2_{\Delta\nu}>\Delta\chi^2)$, "
        r"$\Delta\nu=n_{\rm free}$ (D, $m_c$, [$\ell_0$], [shape]); "
        r"$\sigma=\Phi^{-1}(1-p)$ — same method as `catwise-sigma-null.py`. "
        r"Fixing $\ell_0=15.2^\circ$ removes one free parameter."
    )
    lines.append("")
    lines.append(r"| law | $\ell_0$ mode | free params | $n_{\rm free}$ | $\chi^2$ | $\Delta\chi^2$ | $p$ | $\sigma$ nested | $\sqrt{\Delta\chi^2}$ | $R^2$ |")
    lines.append("|---|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for fp in sorted(fits.values(), key=lambda x: (x["kind"], x["mode"])):
        lines.append(
            f"| {fp['kind']} | {fp['mode']} | {','.join(fp.get('free_params', []))} | "
            f"{fp['n_free']} | {fp['chi2']:.3f} | {fp['Delta_chi2_vs_null']:.3f} | "
            f"{fp['p_value_nested']:.3e} | {fp['sigma_nested_one_sided']:.2f} | "
            f"{fp['sigma_sqrt_Delta_chi2']:.2f} | {fp['R2_unweighted']:.4f} |"
        )
    lines.append("")
    lines.append("## Comparison priors")
    lines.append("")
    lines.append(
        rf"| model | $R^2$ | $\chi^2$ | $\ell_0$ | note |"
    )
    lines.append("|---|---:|---:|---:|---|")
    lines.append(
        f"| ptshell WLS | {PTSHELL_WLS['R2_unweighted']:.3f} | "
        f"{PTSHELL_WLS['chi2']:.1f} | {PTSHELL_WLS['l0']:.1f} | geometry ptshell |"
    )
    lines.append(
        f"| CatWISE-$L$ exp | {CATWISE_EXP_D14['R2_unweighted']:.3f} | "
        f"{CATWISE_EXP_D14['chi2']:.1f} | {CATWISE_EXP_D14['l0']:.1f} | "
        r"$\delta$ lock, $\eta=1$ |"
    )
    lines.append(
        f"| CatWISE-η edge | {CATWISE_ETA_EDGE['R2_unweighted']:.3f} | "
        f"{CATWISE_ETA_EDGE['chi2']:.1f} | {CATWISE_ETA_EDGE['l0']:.1f} | "
        r"η≈0.05 cliff — Chris dislikes |"
    )
    lines.append(
        f"| free-exp (no δ lock) | {FREE_EXP_WLS['R2_unweighted']:.3f} | "
        f"{FREE_EXP_WLS['chi2']:.1f} | — | over-produces contrast |"
    )
    lines.append("")
    lines.append("## Recommendation")
    lines.append("")
    # pick best gentle under fixed15.2, else free
    gentle_fixed = [
        fp for fp in fits.values()
        if fp.get("mode") == "fixed15.2" and fp.get("cliff_index", 99) < 1.5
    ]
    gentle_fixed.sort(key=lambda x: x["chi2"])
    gentle_free = [
        fp for fp in fits.values()
        if fp.get("mode") == "free" and fp.get("cliff_index", 99) < 1.5
    ]
    gentle_free.sort(key=lambda x: x["chi2"])
    if gentle_fixed:
        g = gentle_fixed[0]
        lines.append(
            f"Under preferred $\\ell_0={L0_PTSHELL}^\\circ$, the gentlest "
            f"competitive law is **{g['kind']}** "
            f"($R^2={g['R2_unweighted']:.4f}$, $\\chi^2={g['chi2']:.3f}$, "
            f"cliff={g['cliff_index']:.3g}). "
            + (
                f"Beats CatWISE-exp $R^2\\approx 0.72$ with $\\ell_0\\le 16^\\circ$."
                if g["R2_unweighted"] > 0.72
                else "Does **not** beat CatWISE-exp 0.72 at fixed $\\ell_0=15.2^\\circ$."
            )
        )
    lines.append("")
    if gentle_free:
        g = gentle_free[0]
        lines.append(
            f"With free $\\ell_0$, best gentle law: **{g['kind']}** "
            f"($R^2={g['R2_unweighted']:.4f}$, $\\ell_0={g['l0']:.1f}^\\circ$)."
        )
    lines.append("")
    lines.append(
        "η-edge (η≈0.05) is excluded from the winner set: it creates a "
        "precipice near Earth (cliff-then-flat), which is what Chris wants "
        "to avoid. Soft core+tail / cored / gauss / pl are the physically "
        "gentle candidates."
    )
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("| Path | Content |")
    lines.append("|---|---|")
    lines.append(
        "| `catwise-density-levelset.png` / `.csv` / `.md` | "
        r"$\rho(r)$ overlay, all $\delta=3\%$ laws |"
    )
    lines.append(
        "| `catwise-density-levelset-fit.json` | machine-readable fits |"
    )
    lines.append(
        "| `catwise-density-levelset-migkas-*.png` | per-law / overlay Migkas curves |"
    )
    lines.append(
        "| `catwise-density-levelset-migkas.py` | this script |"
    )
    lines.append(
        "| `prof_core.py` (bubble) | soft kind: "
        r"`(1+|x|^α)^{-β}` |"
    )
    lines.append("")

    md = OUT / "catwise-density-levelset.md"
    md.write_text("\n".join(lines))
    print("wrote", md)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=== CatWISE density level-set + Migkas ===", flush=True)
    print(f"d_lock = {D_LOCK:.4f} Gpc; δ = {DELTA_USE}; D_Earth = {D_EARTH}",
          flush=True)

    # Sanity: δ recovery for each law at D=58.2
    print("\n--- δ-lock sanity (D=58.2) ---", flush=True)
    sanity_shapes = {
        "exp": {},
        "pl": {"n": 3.0},
        "cored": {},
        "gauss": {},
        "soft": {"alpha": 1.5, "beta": SOFT_BETA},
        "eta": {"eta": 0.5},
    }
    for kind, sh in sanity_shapes.items():
        c = density_card(kind, D_EARTH, D_LOCK, sh)
        print(
            f"  {kind:6s} L={c['L_Gpc']:.4g}  δ={c['delta_check']:.6f}  "
            f"relerr={c['delta_rel_err']:.3e}  cliff={c['cliff_index']:.3g}  "
            f"ρ0/ρc={c['rho_0_over_rho_crit']:.4g}",
            flush=True,
        )
        assert c["delta_rel_err"] < 0.01, f"δ miss for {kind}"

    # Migkas fits
    laws = ["exp", "pl", "cored", "gauss", "soft"]
    # optional η≥0.25 — run and keep only if improves over exp
    optional_eta = True

    fits = {}
    curves = {}
    for kind in laws:
        for mode in ("free", "fixed15.2"):
            print(f"\n=== Fit {kind} / {mode} ===", flush=True)
            fitp, full = fit_law(kind, D_LOCK, l0_mode=mode)
            key = f"{kind}_{mode}"
            print("BEST", {k: (round(v, 6) if isinstance(v, float) else v)
                           for k, v in fitp.items()
                           if k not in ("fixed_L_rule",)}, flush=True)
            fits[key] = fitp
            curves[key] = full
            plot_migkas(fitp, full, tag=key)

    if optional_eta:
        print("\n=== Optional η≥0.25 (skip edge) ===", flush=True)
        eta_keep = False
        for mode in ("free", "fixed15.2"):
            try:
                fitp, full = fit_law("eta", D_LOCK, l0_mode=mode)
            except Exception as e:
                print(f"  eta/{mode} failed: {e}", flush=True)
                continue
            exp_key = f"exp_{mode}"
            exp_chi2 = fits[exp_key]["chi2"]
            improves = fitp["chi2"] < exp_chi2 - 0.05 and fitp["shape"]["eta"] >= ETA_FLOOR
            print(
                f"  eta/{mode}: χ²={fitp['chi2']:.3f} vs exp {exp_chi2:.3f}; "
                f"η={fitp['shape']['eta']:.3g}; keep={improves}",
                flush=True,
            )
            if improves:
                key = f"eta_{mode}"
                fits[key] = fitp
                curves[key] = full
                plot_migkas(fitp, full, tag=key)
                eta_keep = True
        if not eta_keep:
            print("  (η≥0.25 does not improve enough — omitted from winner set)",
                  flush=True)

    # Overlays per mode
    for mode in ("free", "fixed15.2"):
        subset = {
            k.split("_")[0]: (fits[k], curves[k])
            for k in fits
            if fits[k]["mode"] == mode
        }
        if subset:
            plot_migkas_overlay(subset, mode)

    # Best shapes for density plot (prefer fixed15.2 if available, else free)
    best_shapes = {}
    for kind in laws + (["eta"] if any(k.startswith("eta_") for k in fits) else []):
        for mode in ("fixed15.2", "free"):
            key = f"{kind}_{mode}"
            if key in fits and fits[key].get("shape"):
                best_shapes[kind] = fits[key]["shape"]
                break
        if kind not in best_shapes:
            if kind == "pl":
                best_shapes[kind] = {"n": 3.0}
            elif kind == "soft":
                best_shapes[kind] = {"alpha": 1.5, "beta": SOFT_BETA}
            elif kind == "eta":
                best_shapes[kind] = {"eta": 0.5}
            else:
                best_shapes[kind] = {}

    print("\n=== Density overlay ===", flush=True)
    dens = plot_density_levelset({}, best_shapes=best_shapes)

    write_md(dens, fits)

    # JSON summary
    payload = dict(
        formula=(
            "level set of smooth ρ(r) with δ=[ρ(D−d)−ρ(D+d)]/ρ(D)=0.031; "
            "A/A_all=a_f(ℓ−ℓ₀)/⟨a_f⟩; Ob=0; L locked; free {D,mc,[ℓ₀],[shape]}"
        ),
        delta_target=DELTA_USE,
        D_CatWISE=D_CATWISE,
        D_Earth_Gpc_density_plot=D_EARTH,
        d_lock_Gpc=float(D_LOCK),
        Omega_m=OMEGA_M,
        R_vision_Gpc=float(RV),
        l0_ptshell=L0_PTSHELL,
        chi2_null=CHI2_NULL,
        N_bins=N_BINS,
        soft_beta_fixed=SOFT_BETA,
        eta_floor=ETA_FLOOR,
        density_params=dens,
        migkas_fits=fits,
        comparison=dict(
            ptshell_wls=PTSHELL_WLS,
            catwise_exp_d14=CATWISE_EXP_D14,
            catwise_eta_edge=CATWISE_ETA_EDGE,
            free_exp_wls=FREE_EXP_WLS,
        ),
        prof_core_soft="kind 'soft': (1+|x|^alpha)^(-beta); LB alpha,beta",
        note=(
            "Primary CatWISE lock uses d=vision≈14.3 Gpc. "
            "Density plot uses D=58.2 Gpc and best-fit free-shape params. "
            "ℓ₀ modes: free and fixed 15.2°. "
            "η-edge (η≈0.05) excluded; optional η≥0.25 kept only if improves."
        ),
    )
    # Strip huge arrays if any accidentally present
    jpath = OUT / "catwise-density-levelset-fit.json"
    jpath.write_text(json.dumps(payload, indent=1, default=str))
    print("wrote", jpath, flush=True)

    # Final ranking printout
    print("\n======== RANKING (ℓ₀=15.2° fixed) ========", flush=True)
    rows = [fp for fp in fits.values() if fp["mode"] == "fixed15.2"]
    rows.sort(key=lambda x: x["chi2"])
    for i, fp in enumerate(rows, 1):
        print(
            f"  {i}. {fp['kind']:6s}  R²={fp['R2_unweighted']:.4f}  "
            f"χ²={fp['chi2']:.3f}  ℓ₀={fp['l0']:.2f}  "
            f"cliff={fp['cliff_index']:.3g}  σ={fp['sigma_nested_one_sided']:.2f} (ν={fp['n_free']})  "
            f"beats0.72@ℓ₀≤16={fp['R2_unweighted']>0.72 and fp['l0']<=16}",
            flush=True,
        )
    print("\n======== RANKING (free ℓ₀) ========", flush=True)
    rows = [fp for fp in fits.values() if fp["mode"] == "free"]
    rows.sort(key=lambda x: x["chi2"])
    for i, fp in enumerate(rows, 1):
        print(
            f"  {i}. {fp['kind']:6s}  R²={fp['R2_unweighted']:.4f}  "
            f"χ²={fp['chi2']:.3f}  ℓ₀={fp['l0']:.2f}  "
            f"cliff={fp['cliff_index']:.3g}  σ={fp['sigma_nested_one_sided']:.2f} (ν={fp['n_free']})  "
            f"beats0.72@ℓ₀≤16={fp['R2_unweighted']>0.72 and fp['l0']<=16}",
            flush=True,
        )


if __name__ == "__main__":
    main()

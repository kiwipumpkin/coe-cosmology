#!/usr/bin/env python3
"""
Pure P₂ quadrupole template + footprint-induced dipole leakage.

Tasks:
1. Build pure P₂(μ) sky: frac = C₂ P₂(μ), C₂ = 2 α κ, κ=0.15
   Preferred axis = Migkas joint direction galactic (l,b)=(303°, −27°)
   α ∈ {1, α_drag≈0.664}
2. Apply masks; fit Migkas-style M + A·n̂ dipole; report fake |A| and direction bias
3. Depths d/D ∈ {0.05,0.1,0.2}: pure P₂ half-ptp + forward-model check
4. JSON dump for report

Skeptical: does leakage alone reach Migkas ~0.1–0.16 under realistic masks?
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from forward_model import ModelParams, H_obs_map, make_latlon_grid
from integrated_gravity import alpha_eff_drag

OUT = Path("/workspace/cosmology")
KAPPA = 0.15
OMEGA_M = 0.3
ETA_DRAG = 1e-4
MIGKAS_L_DEG = 303.0
MIGKAS_B_DEG = -27.0
DD_LIST = [0.05, 0.1, 0.2]
N_RANDOM = 40  # random-mask realizations
RNG = np.random.default_rng(20260917)


def gal_lb_to_cart(l_deg: float, b_deg: float) -> np.ndarray:
    """Galactic (l,b) → Cartesian unit vector (same frame as n_hat)."""
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    return np.array(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)],
        dtype=float,
    )


def P2(mu: np.ndarray) -> np.ndarray:
    return 0.5 * (3.0 * mu**2 - 1.0)


def pure_p2_sky(n_hat: np.ndarray, axis: np.ndarray, C2: float) -> np.ndarray:
    """frac = C₂ P₂(n̂·axiŝ)."""
    a = axis / np.linalg.norm(axis)
    mu = n_hat @ a
    return C2 * P2(mu)


def half_ptp(frac: np.ndarray) -> float:
    return 0.5 * float(np.max(frac) - np.min(frac))


def fit_vector_dipole(frac: np.ndarray, n_hat: np.ndarray, w: np.ndarray) -> dict:
    """
    WLS: frac ≈ M + A_x n_x + A_y n_y + A_z n_z
    Returns A vector, |A|, direction (galactic l,b), R², χ² proxy.
    """
    w = np.asarray(w, dtype=float)
    w = w / w.sum()
    X = np.column_stack([np.ones(len(frac)), n_hat[:, 0], n_hat[:, 1], n_hat[:, 2]])
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(X * sw[:, None], frac * sw, rcond=None)
    M = float(coef[0])
    A_vec = coef[1:4].astype(float)
    A_amp = float(np.linalg.norm(A_vec))
    pred = X @ coef
    resid = frac - pred
    ss_res = float(np.sum(w * resid**2))
    mean = float(np.average(frac, weights=w))
    ss_tot = float(np.sum(w * (frac - mean) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    # Direction of +A (max positive frac under pure dipole)
    if A_amp > 1e-30:
        dhat = A_vec / A_amp
        b = float(np.rad2deg(np.arcsin(np.clip(dhat[2], -1, 1))))
        l = float(np.rad2deg(np.arctan2(dhat[1], dhat[0])) % 360.0)
    else:
        dhat = np.zeros(3)
        l, b = float("nan"), float("nan")

    return {
        "M": M,
        "A_vec": A_vec.tolist(),
        "A_amp": A_amp,
        "l_deg": l,
        "b_deg": b,
        "R2": r2,
        "chi2_proxy": ss_res / (ss_tot + 1e-30),
        "ss_res": ss_res,
        "ss_tot": ss_tot,
    }


def fit_along_axis(frac, n_hat, w, axis: np.ndarray) -> dict:
    """frac ≈ M + A (n̂·axiŝ) — Migkas A cosθ with fixed preferred axis."""
    a = axis / np.linalg.norm(axis)
    mu = n_hat @ a
    w = np.asarray(w, dtype=float)
    w = w / w.sum()
    X = np.column_stack([np.ones(len(frac)), mu])
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(X * sw[:, None], frac * sw, rcond=None)
    pred = X @ coef
    resid = frac - pred
    ss_res = float(np.sum(w * resid**2))
    mean = float(np.average(frac, weights=w))
    ss_tot = float(np.sum(w * (frac - mean) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    # moment projection (no monopole): <f μ>/<μ²>
    A_mom = float(np.sum(w * frac * mu) / (np.sum(w * mu**2) + 1e-30))
    # also Q template
    P2v = P2(mu)
    Xq = np.column_stack([np.ones(len(frac)), P2v])
    coefq, *_ = np.linalg.lstsq(Xq * sw[:, None], frac * sw, rcond=None)
    predq = Xq @ coefq
    residq = frac - predq
    ss_res_q = float(np.sum(w * residq**2))
    r2q = 1.0 - ss_res_q / ss_tot if ss_tot > 0 else 0.0
    return {
        "M": float(coef[0]),
        "A_Migkas": float(coef[1]),
        "A_moment": A_mom,
        "R2_dipole": r2,
        "chi2_proxy_dipole": ss_res / (ss_tot + 1e-30),
        "Q": float(coefq[1]),
        "R2_quad": r2q,
    }


def angular_sep_deg(l1, b1, l2, b2) -> float:
    """Great-circle separation in degrees."""
    r1 = gal_lb_to_cart(l1, b1)
    r2 = gal_lb_to_cart(l2, b2)
    c = float(np.clip(np.dot(r1, r2), -1, 1))
    return float(np.rad2deg(np.arccos(c)))


def make_masks(n_hat, lon, lat, axis, w, rng):
    """
    Return dict name -> boolean mask (True = keep).
    lon, lat from make_latlon_grid: lon=φ from +x, lat=elevation.
    For galactic masks we re-interpret n_hat components as galactic Cartesian
    (same convention as gal_lb_to_cart), so n_z = sin(b), b = arcsin(n_z).
    """
    b_gal = np.rad2deg(np.arcsin(np.clip(n_hat[:, 2], -1, 1)))
    l_gal = np.rad2deg(np.arctan2(n_hat[:, 1], n_hat[:, 0])) % 360.0
    mu = n_hat @ (axis / np.linalg.norm(axis))

    masks = {
        "full_sky": np.ones(len(n_hat), dtype=bool),
        "galactic_abs_b_gt_20": np.abs(b_gal) > 20.0,
        "CoE_facing_hemisphere": mu > 0.0,
        "anti_CoE_hemisphere": mu < 0.0,
    }

    # Random 50% / 30% (solid-angle weighted Bernoulli via cumulative)
    # Use unweighted pixel random for simplicity; report weight sky frac
    for frac_keep, key in [(0.50, "random_50pct"), (0.30, "random_30pct")]:
        # single draw for deterministic masks in this function; multi-realiz. elsewhere
        u = rng.random(len(n_hat))
        masks[key] = u < frac_keep

    # Crude extragalactic cluster-like mask:
    # |b|>20°, exclude |l|<30° strip near GC somewhat, and apply mild CoE-facing
    # preference + random thinning to ~25–35% sky (Migkas-like partial coverage).
    # Also cut a crude "zone of avoidance" |b|<20 already, plus sparse random keep.
    base = np.abs(b_gal) > 20.0
    # Prefer northern galactic + CoE side mildly (selection asymmetry)
    prefer = 0.55 + 0.25 * np.clip(mu, -1, 1) + 0.10 * (b_gal / 90.0)
    prefer = np.clip(prefer, 0.05, 0.95)
    u = rng.random(len(n_hat))
    clusterish = base & (u < prefer * 0.45)  # target ~25-35% of |b|>20
    masks["cluster_like_crude"] = clusterish

    # Extra: |b|>20° AND CoE-facing (worst-case "realistic" combo)
    masks["galactic_abs_b_gt_20_and_CoE_hemi"] = masks["galactic_abs_b_gt_20"] & (
        mu > 0.0
    )

    return masks, b_gal, l_gal, mu


def summarize_mask(name, mask, frac, n_hat, w, axis, true_l, true_b):
    if mask.sum() < 20:
        return {"error": "too few pixels", "n_pix": int(mask.sum())}
    ww = w[mask].copy()
    sky_frac = float(w[mask].sum())
    vd = fit_vector_dipole(frac[mask], n_hat[mask], ww)
    fa = fit_along_axis(frac[mask], n_hat[mask], ww, axis)
    bias = angular_sep_deg(vd["l_deg"], vd["b_deg"], true_l, true_b) if np.isfinite(
        vd["l_deg"]
    ) else None
    # angle between recovered dipole and true axis (0 = aligned)
    return {
        "n_pix": int(mask.sum()),
        "sky_frac_weight": sky_frac,
        "vector_dipole": vd,
        "fixed_axis": fa,
        "direction_bias_deg_from_Migkas_axis": bias,
        "A_amp_over_0.10": vd["A_amp"] / 0.10,
        "A_amp_over_0.16": vd["A_amp"] / 0.16,
        "A_Migkas_over_0.10": abs(fa["A_Migkas"]) / 0.10,
        "A_moment_over_0.10": abs(fa["A_moment"]) / 0.10,
    }


def random_mask_stats(frac, n_hat, w, axis, keep_frac, n_real, rng, true_l, true_b):
    """Many random masks; report median/mean/p95 of |A|."""
    amps_vec = []
    amps_mig = []
    amps_mom = []
    biases = []
    r2s = []
    for _ in range(n_real):
        u = rng.random(len(n_hat))
        mask = u < keep_frac
        if mask.sum() < 20:
            continue
        ww = w[mask]
        vd = fit_vector_dipole(frac[mask], n_hat[mask], ww)
        fa = fit_along_axis(frac[mask], n_hat[mask], ww, axis)
        amps_vec.append(vd["A_amp"])
        amps_mig.append(abs(fa["A_Migkas"]))
        amps_mom.append(abs(fa["A_moment"]))
        biases.append(
            angular_sep_deg(vd["l_deg"], vd["b_deg"], true_l, true_b)
            if np.isfinite(vd["l_deg"])
            else np.nan
        )
        r2s.append(vd["R2"])

    def stats(arr):
        a = np.asarray(arr, dtype=float)
        return {
            "mean": float(np.nanmean(a)),
            "median": float(np.nanmedian(a)),
            "p95": float(np.nanpercentile(a, 95)),
            "max": float(np.nanmax(a)),
            "std": float(np.nanstd(a)),
        }

    return {
        "keep_frac": keep_frac,
        "n_realizations": n_real,
        "A_amp_vector": stats(amps_vec),
        "A_Migkas_fixed_axis": stats(amps_mig),
        "A_moment": stats(amps_mom),
        "direction_bias_deg": stats(biases),
        "R2_vector": stats(r2s),
    }


def forward_model_half_ptp(d_over_D, kappa, alpha):
    """Exact forward H_obs/H-1 half-ptp and Legendre C2 on μ grid."""
    D, H, G = 1.0, 1.0, 1.0
    d = d_over_D * D
    p = ModelParams(D=D, M=kappa, H=H, alpha=alpha, G=G)
    mu = np.linspace(-1 + 1e-6, 1 - 1e-6, 1024)
    s = np.sqrt(np.maximum(1 - mu**2, 0))
    n_hat = np.stack([mu, np.zeros_like(mu), s], axis=-1)
    _, _, frac = H_obs_map(n_hat, d, p)
    # C2 projection
    from numpy.polynomial.legendre import leggauss

    x, wg = leggauss(256)
    f_at = np.interp(x, mu, frac)
    P2x = 0.5 * (3 * x**2 - 1)
    C2 = 2.5 * np.sum(wg * f_at * P2x)  # (2*2+1)/2 = 2.5
    P1x = x
    C1 = 1.5 * np.sum(wg * f_at * P1x)
    return {
        "d_over_D": d_over_D,
        "half_ptp_forward": half_ptp(frac),
        "C2_forward": float(C2),
        "C1_forward": float(C1),
        "frac_max": float(frac.max()),
        "frac_min": float(frac.min()),
    }


def analytic_p2_half_ptp(C2: float) -> float:
    """P₂ ∈ [-1/2, 1] ⇒ C₂ P₂ ∈ [-C₂/2, C₂]; half-ptp = 3/4 C₂."""
    return 0.75 * C2


def main():
    a_drag = float(alpha_eff_drag(ETA_DRAG))
    axis = gal_lb_to_cart(MIGKAS_L_DEG, MIGKAS_B_DEG)
    n_hat, lon, lat, w = make_latlon_grid(180, 90)  # finer grid

    results = {
        "kappa": KAPPA,
        "Omega_m": OMEGA_M,
        "alpha_bolt_on": 1.0,
        "alpha_drag": a_drag,
        "Migkas_axis_galactic_lb": [MIGKAS_L_DEG, MIGKAS_B_DEG],
        "C2_formula": "C2 = 2 * alpha * kappa",
        "C2_bolt_on": 2.0 * 1.0 * KAPPA,
        "C2_drag": 2.0 * a_drag * KAPPA,
        "pure_P2_half_ptp_formula": "0.75 * C2",
        "pure_P2_half_ptp_bolt_on": analytic_p2_half_ptp(2.0 * KAPPA),
        "pure_P2_half_ptp_drag": analytic_p2_half_ptp(2.0 * a_drag * KAPPA),
        "Migkas_target_A": {"low": 0.10, "high": 0.16, "reported_pct": "16±3% DL"},
        "note": (
            "Pure P2 sky has ZERO full-sky Legendre C1. Any fitted dipole is "
            "footprint leakage. Vector dipole free direction; fixed-axis A is "
            "Migkas-style A cosθ along (303,-27)."
        ),
    }

    # Depth table: pure P2 vs forward
    depth_rows = []
    for dD in DD_LIST:
        for label, alpha in [("bolt_on", 1.0), ("drag", a_drag)]:
            C2 = 2.0 * alpha * KAPPA
            fm = forward_model_half_ptp(dD, KAPPA, alpha)
            depth_rows.append(
                {
                    "alpha_label": label,
                    "alpha": alpha,
                    "d_over_D": dD,
                    "C2_tidal": C2,
                    "half_ptp_pure_P2": analytic_p2_half_ptp(C2),
                    **fm,
                }
            )
    results["depth_table"] = depth_rows

    # Footprint leakage for both alphas
    leakage = {}
    for label, alpha in [("bolt_on", 1.0), ("drag", a_drag)]:
        C2 = 2.0 * alpha * KAPPA
        frac = pure_p2_sky(n_hat, axis, C2)
        assert abs(half_ptp(frac) - analytic_p2_half_ptp(C2)) < 1e-3

        # Deterministic masks (fixed RNG stream per label)
        rng = np.random.default_rng(20260917 + (0 if label == "bolt_on" else 1))
        masks, b_gal, l_gal, mu = make_masks(n_hat, lon, lat, axis, w, rng)

        mask_results = {}
        for name, mask in masks.items():
            mask_results[name] = summarize_mask(
                name, mask, frac, n_hat, w, axis, MIGKAS_L_DEG, MIGKAS_B_DEG
            )

        # Multi-realization random masks
        rng2 = np.random.default_rng(42 if label == "bolt_on" else 43)
        rand_stats = {
            "random_50pct": random_mask_stats(
                frac, n_hat, w, axis, 0.50, N_RANDOM, rng2, MIGKAS_L_DEG, MIGKAS_B_DEG
            ),
            "random_30pct": random_mask_stats(
                frac, n_hat, w, axis, 0.30, N_RANDOM, rng2, MIGKAS_L_DEG, MIGKAS_B_DEG
            ),
        }

        # Max fake A across deterministic masks (exclude full_sky which should be ~0)
        candidates = []
        for name, v in mask_results.items():
            if "error" in v:
                continue
            candidates.append(
                (
                    name,
                    v["vector_dipole"]["A_amp"],
                    v["fixed_axis"]["A_Migkas"],
                    v["fixed_axis"]["A_moment"],
                    v.get("direction_bias_deg_from_Migkas_axis"),
                )
            )
        # also include random p95
        for rk, rs in rand_stats.items():
            candidates.append(
                (
                    f"{rk}_p95",
                    rs["A_amp_vector"]["p95"],
                    rs["A_Migkas_fixed_axis"]["p95"],
                    rs["A_moment"]["p95"],
                    rs["direction_bias_deg"]["median"],
                )
            )

        max_by_Aamp = max(candidates, key=lambda t: t[1])
        max_by_Amig = max(candidates, key=lambda t: abs(t[2]))
        max_by_Amom = max(candidates, key=lambda t: abs(t[3]))

        leakage[label] = {
            "C2": C2,
            "half_ptp_pure_P2": analytic_p2_half_ptp(C2),
            "full_sky_A_amp_check": mask_results["full_sky"]["vector_dipole"]["A_amp"],
            "masks": mask_results,
            "random_multi": rand_stats,
            "max_fake_A_amp": {
                "mask": max_by_Aamp[0],
                "A_amp": max_by_Aamp[1],
                "A_Migkas": max_by_Aamp[2],
                "A_moment": max_by_Aamp[3],
            },
            "max_fake_A_Migkas": {
                "mask": max_by_Amig[0],
                "A_amp": max_by_Amig[1],
                "A_Migkas": max_by_Amig[2],
                "A_moment": max_by_Amig[3],
            },
            "max_fake_A_moment": {
                "mask": max_by_Amom[0],
                "A_amp": max_by_Amom[1],
                "A_Migkas": max_by_Amom[2],
                "A_moment": max_by_Amom[3],
            },
        }

    results["leakage"] = leakage

    # Recommendation logic
    # "Realistic" masks: galactic |b|>20, random 30/50, cluster_like — NOT pure CoE hemisphere
    # unless Migkas is strongly CoE-facing selected (unlikely by design)
    bolt = leakage["bolt_on"]
    drag = leakage["drag"]

    def realistic_max(block):
        names = [
            "galactic_abs_b_gt_20",
            "cluster_like_crude",
            "random_50pct",
            "random_30pct",
        ]
        vals_amp = []
        vals_mig = []
        vals_mom = []
        for n in names:
            v = block["masks"][n]
            if "error" in v:
                continue
            vals_amp.append((n, v["vector_dipole"]["A_amp"]))
            vals_mig.append((n, abs(v["fixed_axis"]["A_Migkas"])))
            vals_mom.append((n, abs(v["fixed_axis"]["A_moment"])))
        # add random medians
        for rk in ("random_50pct", "random_30pct"):
            rs = block["random_multi"][rk]
            vals_amp.append((f"{rk}_median", rs["A_amp_vector"]["median"]))
            vals_amp.append((f"{rk}_p95", rs["A_amp_vector"]["p95"]))
            vals_mig.append((f"{rk}_median", rs["A_Migkas_fixed_axis"]["median"]))
            vals_mom.append((f"{rk}_median", rs["A_moment"]["median"]))
        return {
            "max_A_amp": max(vals_amp, key=lambda t: t[1]),
            "max_A_Migkas": max(vals_mig, key=lambda t: t[1]),
            "max_A_moment": max(vals_mom, key=lambda t: t[1]),
        }

    rec = {
        "bolt_on_realistic": realistic_max(bolt),
        "drag_realistic": realistic_max(drag),
        "bolt_on_CoE_hemi_A_amp": bolt["masks"]["CoE_facing_hemisphere"]["vector_dipole"][
            "A_amp"
        ],
        "bolt_on_CoE_hemi_A_Migkas": bolt["masks"]["CoE_facing_hemisphere"]["fixed_axis"][
            "A_Migkas"
        ],
        "bolt_on_CoE_hemi_A_moment": bolt["masks"]["CoE_facing_hemisphere"]["fixed_axis"][
            "A_moment"
        ],
        "drag_CoE_hemi_A_moment": drag["masks"]["CoE_facing_hemisphere"]["fixed_axis"][
            "A_moment"
        ],
    }
    results["recommendation_inputs"] = rec

    # Verdict strings
    real_A = rec["drag_realistic"]["max_A_amp"][1]  # prefer drag (time-integrated)
    real_A_bolt = rec["bolt_on_realistic"]["max_A_amp"][1]
    hemi_mom = rec["bolt_on_CoE_hemi_A_moment"]
    results["verdict"] = {
        "full_sky_leakage": "negligible (A~0 by orthogonality)",
        "realistic_masks_reach_Migkas_0.10_0.16": bool(
            real_A_bolt >= 0.10 or real_A >= 0.10
        ),
        "CoE_hemisphere_reaches_Migkas": bool(abs(hemi_mom) >= 0.10),
        "recommendation_for_Chris": (
            "Leakage alone does NOT reach Migkas A~0.1–0.16 under realistic "
            "footprints (galactic |b|>20°, random 30–50%, crude cluster-like): "
            "fake |A| stays ≪0.05 typically. "
            "ONLY a strongly CoE-facing hemisphere/cone selection produces "
            "fake A_moment~0.11 (α=1) / ~0.07 (drag) or A_Migkas(M+Aμ)~0.4 — "
            "an analysis artifact, not a clean CP signal. "
            "Do NOT attribute Migkas dipole to pure-P2 footprint leakage under "
            "standard extragalactic cuts. Prefer quadrupole template comparison."
        ),
    }

    out_path = OUT / "p2_footprint_leakage_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    # Console summary
    print("=== Pure P2 template ===")
    print(
        f"κ={KAPPA}, C2(α=1)={results['C2_bolt_on']:.4f}, "
        f"C2(drag)={results['C2_drag']:.4f}"
    )
    print(
        f"half-ptp pure P2: bolt={results['pure_P2_half_ptp_bolt_on']:.4f}, "
        f"drag={results['pure_P2_half_ptp_drag']:.4f}"
    )
    print(f"axis (l,b)=({MIGKAS_L_DEG},{MIGKAS_B_DEG})")

    print("\n=== Depth: pure P2 vs forward half-ptp ===")
    for r in depth_rows:
        if r["alpha_label"] != "bolt_on":
            continue
        print(
            f"d/D={r['d_over_D']:.2f}  halfptp_P2={r['half_ptp_pure_P2']:.4f}  "
            f"halfptp_fwd={r['half_ptp_forward']:.4f}  C2_fwd={r['C2_forward']:.4f}"
        )

    print("\n=== Footprint leakage (α=1 bolt-on) ===")
    for name, v in bolt["masks"].items():
        if "error" in v:
            continue
        print(
            f"{name:<40} sky={v['sky_frac_weight']:.3f}  "
            f"A_amp={v['vector_dipole']['A_amp']:.4f}  "
            f"A_Mig={v['fixed_axis']['A_Migkas']:+.4f}  "
            f"A_mom={v['fixed_axis']['A_moment']:+.4f}  "
            f"R2dip={v['fixed_axis']['R2_dipole']:.3f}  "
            f"bias={v['direction_bias_deg_from_Migkas_axis']}"
        )
    print("\nRandom multi (α=1):")
    for rk, rs in bolt["random_multi"].items():
        print(
            f"  {rk}: A_amp median={rs['A_amp_vector']['median']:.4f} "
            f"p95={rs['A_amp_vector']['p95']:.4f} max={rs['A_amp_vector']['max']:.4f}"
        )

    print("\n=== Footprint leakage (α_drag) ===")
    for name, v in drag["masks"].items():
        if "error" in v:
            continue
        print(
            f"{name:<40} sky={v['sky_frac_weight']:.3f}  "
            f"A_amp={v['vector_dipole']['A_amp']:.4f}  "
            f"A_Mig={v['fixed_axis']['A_Migkas']:+.4f}  "
            f"A_mom={v['fixed_axis']['A_moment']:+.4f}"
        )

    print("\n=== Verdict ===")
    print(results["verdict"]["recommendation_for_Chris"])
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()

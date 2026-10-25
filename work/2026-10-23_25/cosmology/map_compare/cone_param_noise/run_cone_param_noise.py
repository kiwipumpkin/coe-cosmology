#!/usr/bin/env python3
"""
A/A_i maps with parameter variance INSIDE the calculation (not pixel patchwork).

TWO families (both delivered):

1) CONE-NATIVE (primary physics product)
   For each θ=75° cone center: sample N_samp sightlines uniformly in solid angle
   inside the cone; each object gets optional radial velocity noise and/or
   distance noise; cone A = mean of A_raw over objects; then mean-renorm.

2) ANGLE / spreadsheet-style companion
   One evaluation per cone-center direction n̂ at observed distance; add
   random to v and/or d; then mean-renorm.

Each family: baseline | vel-only | dist-only | vel+dist combined
(+ louder vel/dist companions). Noise never applied by mixing finished pixels.
"""
from __future__ import annotations

import json
import sys
import time
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

ROOT = Path("/workspace/cosmology")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "density_local_match"))
sys.path.insert(0, str(ROOT / "map_compare"))

from run_density_local_match import (  # noqa: E402
    params_twoscale_local_match,
    peculiar_u,
    D,
    H,
    G,
)
from run_levelset_maps import (  # noqa: E402
    lon_disp_from_ell,
    celestial_equator_lb,
    style,
    VMIN,
    VMAX,
    rotation_matrix_x_to_sky,
    physics_sky_maps,
    build_cone_geometry,
    cone_average,
    plot_card as _plot_card_ref,
    CARD_PX,
)

OUT = Path(__file__).resolve().parent
CONE_REF = ROOT / "map_compare" / "cone75_maps.npz"

THETA_DEG = 75.0
N_SAMP = 64
MIGKAS_L, MIGKAS_B = 280.0, -20.0
BASE = dict(x=5.9, f_u=0.94, alpha=1.0, eps=0.1)

# Seeds
SEED_SIGHTLINES = 42  # shared cone-native sightline draws
SEED_CN_VEL = 43
SEED_CN_DIST = 44
SEED_CN_BOTH_V = 45
SEED_CN_BOTH_D = 46
SEED_ANG_VEL = 143
SEED_ANG_DIST = 144
SEED_ANG_BOTH_V = 145
SEED_ANG_BOTH_D = 146

SIGMA_U_FRAC_PRIMARY = 0.05
SIGMA_U_FRAC_LOUD = 0.10
SIGMA_D_REL_PRIMARY = 0.15
SIGMA_D_REL_LOUD = 0.25
D_CLIP_FRAC = 0.02


def now_local() -> str:
    return datetime.now(ZoneInfo("America/Chicago")).strftime(
        "%Y-%m-%d %I:%M %p %Z (UTC%z)"
    )


def sky_unit_vectors(L_centers, B_centers):
    Lc, Bc = np.meshgrid(L_centers, B_centers)
    l = np.deg2rad(Lc)
    b = np.deg2rad(Bc)
    n = np.stack(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], axis=-1
    )
    return n


def to_model_frame(n_sky, R):
    shape = n_sky.shape
    flat = n_sky.reshape(-1, 3)
    return (R.T @ flat.T).T.reshape(shape)


def orthonormal_basis(c):
    c = np.asarray(c, dtype=float)
    c = c / np.linalg.norm(c)
    a = np.array([0.0, 0.0, 1.0]) if abs(c[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(a, c)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(c, e1)
    return e1, e2, c


def sample_in_cone(c_hat, theta_deg, n_samp, rng):
    """Uniform in solid angle within theta of c_hat. Returns (n_samp, 3)."""
    e1, e2, ez = orthonormal_basis(c_hat)
    cos_max = np.cos(np.deg2rad(theta_deg))
    mu = rng.uniform(cos_max, 1.0, size=n_samp)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n_samp)
    st = np.sqrt(np.maximum(1.0 - mu * mu, 0.0))
    n = (
        (st * np.cos(phi))[:, None] * e1
        + (st * np.sin(phi))[:, None] * e2
        + mu[:, None] * ez
    )
    n /= np.linalg.norm(n, axis=-1, keepdims=True)
    return n


def compute_A_raw(n_hat_model, pmeta, alpha, d, xi_vel=None, sigma_u=0.0):
    """
    n_hat_model: (N, 3) model-frame sightlines
    d: scalar or (N,)
    xi_vel: optional (N,) N(0,1); adds ξ σ_u r̂ to u_r
    Returns A_inv_raw, A_hr_raw each (N,)
    """
    n = np.asarray(n_hat_model, dtype=float)
    N = n.shape[0]
    d_arr = np.broadcast_to(np.asarray(d, dtype=float), (N,)).copy()

    D_vec = np.array([D, 0.0, 0.0])
    r_vec = D_vec + d_arr[:, None] * n
    u_r = peculiar_u(r_vec, pmeta, alpha=alpha)
    u_D = peculiar_u(D_vec.reshape(1, 3), pmeta, alpha=alpha)[0]

    if xi_vel is not None and sigma_u > 0.0:
        r_norm = np.linalg.norm(r_vec, axis=-1, keepdims=True)
        r_hat = r_vec / np.maximum(r_norm, 1e-30)
        u_r = u_r + (np.asarray(xi_vel) * sigma_u)[:, None] * r_hat

    delta_v = np.sum((u_r - u_D) * n, axis=-1)
    H_obs = H + delta_v / d_arr
    H_safe = np.maximum(H_obs, 1e-6 * H)
    A_hr = H_safe / H
    A_inv = 1.0 / (H_safe / H) ** 2
    return A_inv, A_hr


def mean_renorm(Z):
    m = float(np.nanmean(Z))
    return Z / m if m != 0 else Z


def noisy_d(d0, eta, sigma_d_rel):
    d = d0 * (1.0 + eta * sigma_d_rel)
    return np.maximum(d, D_CLIP_FRAC * d0)


def plot_card(L_c, B_c, Z, title, path, show_marks=False):
    _plot_card_ref(L_c, B_c, Z, title, path, show_marks=show_marks)


def prep_mesh(L_c, B_c, Z):
    Lon_c, Lat_c = np.meshgrid(L_c, B_c)
    lon_disp = np.deg2rad(lon_disp_from_ell(Lon_c))
    order = np.argsort(lon_disp[0])
    lon_s = lon_disp[:, order]
    Z_s = Z[:, order]
    lon_centers = lon_s[0]
    mid = 0.5 * (lon_centers[:-1] + lon_centers[1:])
    lon_edges = np.concatenate(
        [
            [lon_centers[0] - (mid[0] - lon_centers[0])],
            mid,
            [lon_centers[-1] + (lon_centers[-1] - mid[-1])],
        ]
    )
    dB = float(np.median(np.diff(B_c))) if len(B_c) > 1 else 5.0
    B_edge = np.concatenate([B_c - dB / 2, [B_c[-1] + dB / 2]])
    Lon_e, Lat_e = np.meshgrid(lon_edges, np.deg2rad(B_edge))
    return Lon_e, Lat_e, Z_s


def draw_equator(ax):
    ell_eq, b_eq = celestial_equator_lb()
    if ell_eq is None:
        return
    ld = np.deg2rad(lon_disp_from_ell(ell_eq))
    bd = np.deg2rad(b_eq)
    jumps = np.where(np.abs(np.diff(ld)) > np.pi / 2)[0]
    starts = np.concatenate([[0], jumps + 1])
    ends = np.concatenate([jumps + 1, [len(ld)]])
    for s, e in zip(starts, ends):
        if e - s < 2:
            continue
        ax.plot(ld[s:e], bd[s:e], "k-", lw=0.8, alpha=0.75, zorder=5)


def plot_strip(entries, path, title):
    style()
    n = len(entries)
    fig, axes = plt.subplots(
        1,
        n,
        figsize=(4.2 * n + 1.0, 3.4),
        subplot_kw={"projection": "mollweide"},
    )
    axes = np.atleast_1d(axes).ravel()
    pcm = None
    for ax, ent in zip(axes, entries):
        Lon_e, Lat_e, Z_s = prep_mesh(ent["L"], ent["B"], ent["Z"])
        pcm = ax.pcolormesh(
            Lon_e, Lat_e, Z_s, shading="flat", cmap="magma",
            norm=Normalize(VMIN, VMAX), rasterized=True,
        )
        ax.grid(True, linestyle=":", alpha=0.35, color="white", lw=0.4)
        draw_equator(ax)
        ax.set_xticks(np.deg2rad([-90, 0, 90]))
        ax.set_xticklabels(["90", "0", "270"], fontsize=8)
        ax.set_yticks([])
        ax.set_title(ent["short_title"], fontsize=9, pad=8)
    fig.subplots_adjust(left=0.02, right=0.90, bottom=0.06, top=0.80, wspace=0.08)
    fig.suptitle(title, fontsize=11)
    if pcm is not None:
        cax = fig.add_axes([0.915, 0.15, 0.012, 0.6])
        cbar = fig.colorbar(pcm, cax=cax)
        cbar.set_label(r"$A/A_i$")
        cbar.set_ticks(np.linspace(VMIN, VMAX, 8))
    fig.savefig(path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def plot_compare_grid(rows, path, title):
    """
    rows: list of dicts with key 'label' and 'maps': list of {L,B,Z,short_title}
    Same columns across rows (baseline|vel|dist|both).
    """
    style()
    nrows = len(rows)
    ncols = len(rows[0]["maps"])
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.0 * ncols + 1.0, 2.8 * nrows + 0.8),
        subplot_kw={"projection": "mollweide"},
    )
    if nrows == 1:
        axes = np.array([axes])
    pcm = None
    for i, row in enumerate(rows):
        for j, ent in enumerate(row["maps"]):
            ax = axes[i, j]
            Lon_e, Lat_e, Z_s = prep_mesh(ent["L"], ent["B"], ent["Z"])
            pcm = ax.pcolormesh(
                Lon_e, Lat_e, Z_s, shading="flat", cmap="magma",
                norm=Normalize(VMIN, VMAX), rasterized=True,
            )
            ax.grid(True, linestyle=":", alpha=0.3, color="white", lw=0.35)
            ax.set_xticks(np.deg2rad([-90, 0, 90]))
            ax.set_xticklabels(["90", "0", "270"], fontsize=7)
            ax.set_yticks([])
            if i == 0:
                ax.set_title(ent["short_title"], fontsize=9, pad=6)
            if j == 0:
                ax.text(
                    -0.08, 0.5, row["label"], transform=ax.transAxes,
                    rotation=90, va="center", ha="right", fontsize=10,
                )
    fig.subplots_adjust(
        left=0.06, right=0.90, bottom=0.04, top=0.88, wspace=0.06, hspace=0.18
    )
    fig.suptitle(title, fontsize=12)
    if pcm is not None:
        cax = fig.add_axes([0.915, 0.12, 0.012, 0.7])
        cbar = fig.colorbar(pcm, cax=cax)
        cbar.set_label(r"$A/A_i$")
        cbar.set_ticks(np.linspace(VMIN, VMAX, 8))
    fig.savefig(path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def span_str(Z):
    return f"[{float(np.nanmin(Z)):.4f}, {float(np.nanmax(Z)):.4f}]"


# ---------------------------------------------------------------------------
# Family builders
# ---------------------------------------------------------------------------

def build_angle_family(n_model, pmeta, alpha, d0, sigma_u_05, sigma_u_10):
    """Spreadsheet-style: one eval per cone-center angle."""
    shape = n_model.shape[:2]
    flat = n_model.reshape(-1, 3)
    nflat = flat.shape[0]

    rng_v = np.random.default_rng(SEED_ANG_VEL)
    rng_d = np.random.default_rng(SEED_ANG_DIST)
    rng_bv = np.random.default_rng(SEED_ANG_BOTH_V)
    rng_bd = np.random.default_rng(SEED_ANG_BOTH_D)

    xi05 = rng_v.normal(0.0, 1.0, size=nflat)
    xi10 = rng_v.normal(0.0, 1.0, size=nflat)
    eta15 = rng_d.normal(0.0, 1.0, size=nflat)
    eta25 = rng_d.normal(0.0, 1.0, size=nflat)
    xi_b = rng_bv.normal(0.0, 1.0, size=nflat)
    eta_b = rng_bd.normal(0.0, 1.0, size=nflat)

    out = {}
    noise = {
        "xi_vel_05": xi05.reshape(shape),
        "xi_vel_10": xi10.reshape(shape),
        "eta_dist_15": eta15.reshape(shape),
        "eta_dist_25": eta25.reshape(shape),
        "xi_both_v": xi_b.reshape(shape),
        "eta_both_d": eta_b.reshape(shape),
    }

    def pack(key, A_inv, A_hr):
        Zi = mean_renorm(A_inv.reshape(shape))
        Zh = mean_renorm(A_hr.reshape(shape))
        out[f"{key}_invsq"] = Zi
        out[f"{key}_hratio"] = Zh
        print(f"  [angle] {key}: invsq {span_str(Zi)}", flush=True)
        return Zi

    A_inv, A_hr = compute_A_raw(flat, pmeta, alpha, d0)
    pack("baseline", A_inv, A_hr)

    A_inv, A_hr = compute_A_raw(flat, pmeta, alpha, d0, xi05, sigma_u_05)
    pack("vel_s05", A_inv, A_hr)
    A_inv, A_hr = compute_A_raw(flat, pmeta, alpha, d0, xi10, sigma_u_10)
    pack("vel_s10", A_inv, A_hr)

    A_inv, A_hr = compute_A_raw(flat, pmeta, alpha, noisy_d(d0, eta15, SIGMA_D_REL_PRIMARY))
    pack("dist_s15", A_inv, A_hr)
    A_inv, A_hr = compute_A_raw(flat, pmeta, alpha, noisy_d(d0, eta25, SIGMA_D_REL_LOUD))
    pack("dist_s25", A_inv, A_hr)

    A_inv, A_hr = compute_A_raw(
        flat, pmeta, alpha, noisy_d(d0, eta_b, SIGMA_D_REL_PRIMARY), xi_b, sigma_u_05
    )
    pack("both_s05_s15", A_inv, A_hr)

    return out, noise


def build_conenative_family(c_sky, R, pmeta, alpha, d0, sigma_u_05, sigma_u_10):
    """
    Cone-native: for each cone center, sample N_SAMP sightlines in the cone,
    apply per-object noise, average A_raw → cone value.
    Shared sightline geometry (SEED_SIGHTLINES); independent noise seeds.
    """
    n_cones = c_sky.shape[0]
    # Pre-draw all sightlines in sky frame, then rotate to model
    rng_sl = np.random.default_rng(SEED_SIGHTLINES)
    # (n_cones, N_SAMP, 3) sky
    n_sky = np.empty((n_cones, N_SAMP, 3), dtype=float)
    for i in range(n_cones):
        n_sky[i] = sample_in_cone(c_sky[i], THETA_DEG, N_SAMP, rng_sl)
    n_model = to_model_frame(n_sky, R)  # (n_cones, N_SAMP, 3)

    rng_v = np.random.default_rng(SEED_CN_VEL)
    rng_d = np.random.default_rng(SEED_CN_DIST)
    rng_bv = np.random.default_rng(SEED_CN_BOTH_V)
    rng_bd = np.random.default_rng(SEED_CN_BOTH_D)

    xi05 = rng_v.normal(0.0, 1.0, size=(n_cones, N_SAMP))
    xi10 = rng_v.normal(0.0, 1.0, size=(n_cones, N_SAMP))
    eta15 = rng_d.normal(0.0, 1.0, size=(n_cones, N_SAMP))
    eta25 = rng_d.normal(0.0, 1.0, size=(n_cones, N_SAMP))
    xi_b = rng_bv.normal(0.0, 1.0, size=(n_cones, N_SAMP))
    eta_b = rng_bd.normal(0.0, 1.0, size=(n_cones, N_SAMP))

    flat = n_model.reshape(-1, 3)

    def cone_mean(A_obj):
        """A_obj (n_cones, N_SAMP) → mean over samp, then will renorm later."""
        return A_obj.mean(axis=1)

    out = {}
    noise = {
        "xi_vel_05": xi05,
        "xi_vel_10": xi10,
        "eta_dist_15": eta15,
        "eta_dist_25": eta25,
        "xi_both_v": xi_b,
        "eta_both_d": eta_b,
        "n_samp": N_SAMP,
    }

    def run_variant(key, d_arr, xi, sigma_u):
        # d_arr, xi: None or (n_cones, N_SAMP)
        if d_arr is None:
            d_use = d0
        else:
            d_use = d_arr.reshape(-1)
        xi_use = None if xi is None else xi.reshape(-1)
        A_inv, A_hr = compute_A_raw(flat, pmeta, alpha, d_use, xi_use, sigma_u)
        Ai = cone_mean(A_inv.reshape(n_cones, N_SAMP))
        Ah = cone_mean(A_hr.reshape(n_cones, N_SAMP))
        return Ai, Ah

    def pack(key, Ai, Ah, shape):
        Zi = mean_renorm(Ai.reshape(shape))
        Zh = mean_renorm(Ah.reshape(shape))
        out[f"{key}_invsq"] = Zi
        out[f"{key}_hratio"] = Zh
        print(f"  [conenative] {key}: invsq {span_str(Zi)}", flush=True)
        return Zi

    # shape filled by caller via reshape — we return flat cone vectors;
    # caller reshapes. Actually return (n_cones,) and let main reshape.
    specs = []

    Ai, Ah = run_variant("baseline", None, None, 0.0)
    out["baseline_invsq_flat"] = Ai
    out["baseline_hratio_flat"] = Ah
    specs.append("baseline")

    Ai, Ah = run_variant("vel_s05", None, xi05, sigma_u_05)
    out["vel_s05_invsq_flat"] = Ai
    out["vel_s05_hratio_flat"] = Ah

    Ai, Ah = run_variant("vel_s10", None, xi10, sigma_u_10)
    out["vel_s10_invsq_flat"] = Ai
    out["vel_s10_hratio_flat"] = Ah

    Ai, Ah = run_variant(
        "dist_s15", noisy_d(d0, eta15, SIGMA_D_REL_PRIMARY), None, 0.0
    )
    out["dist_s15_invsq_flat"] = Ai
    out["dist_s15_hratio_flat"] = Ah

    Ai, Ah = run_variant(
        "dist_s25", noisy_d(d0, eta25, SIGMA_D_REL_LOUD), None, 0.0
    )
    out["dist_s25_invsq_flat"] = Ai
    out["dist_s25_hratio_flat"] = Ah

    Ai, Ah = run_variant(
        "both_s05_s15",
        noisy_d(d0, eta_b, SIGMA_D_REL_PRIMARY),
        xi_b,
        sigma_u_05,
    )
    out["both_s05_s15_invsq_flat"] = Ai
    out["both_s05_s15_hratio_flat"] = Ah

    return out, noise, n_sky


def finalize_cn(cn_raw, shape):
    """Mean-renorm flat cone arrays into shaped maps."""
    keys = ["baseline", "vel_s05", "vel_s10", "dist_s15", "dist_s25", "both_s05_s15"]
    out = {}
    for k in keys:
        Ai = cn_raw[f"{k}_invsq_flat"]
        Ah = cn_raw[f"{k}_hratio_flat"]
        Zi = mean_renorm(Ai.reshape(shape))
        Zh = mean_renorm(Ah.reshape(shape))
        out[f"{k}_invsq"] = Zi
        out[f"{k}_hratio"] = Zh
        print(f"  [conenative] {k}: invsq {span_str(Zi)}", flush=True)
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t_all = time.time()
    print(f"[{now_local()}] cone_param_noise (BOTH families) → {OUT}", flush=True)

    ref = np.load(CONE_REF)
    L_c = np.asarray(ref["L_centers"], dtype=float)
    B_c = np.asarray(ref["B_centers"], dtype=float)
    nB, nL = len(B_c), len(L_c)
    shape = (nB, nL)
    print(f"grid: L={nL} B={nB}  N_samp={N_SAMP}  θ={THETA_DEG}°", flush=True)

    pmeta = params_twoscale_local_match(BASE["x"], BASE["f_u"])
    alpha = BASE["alpha"]
    d0 = BASE["eps"] * D
    R = rotation_matrix_x_to_sky(MIGKAS_L, MIGKAS_B)
    sigma_u_05 = SIGMA_U_FRAC_PRIMARY * H * d0
    sigma_u_10 = SIGMA_U_FRAC_LOUD * H * d0

    n_sky_grid = sky_unit_vectors(L_c, B_c)
    n_model_grid = to_model_frame(n_sky_grid, R)
    Lc, Bc = np.meshgrid(L_c, B_c)
    l = np.deg2rad(Lc.ravel())
    b = np.deg2rad(Bc.ravel())
    c_sky = np.stack(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], axis=-1
    )

    # ----- Family 1: cone-native -----
    print("\n=== Family 1: CONE-NATIVE (multi-sightline) ===", flush=True)
    t0 = time.time()
    cn_raw, cn_noise, _ = build_conenative_family(
        c_sky, R, pmeta, alpha, d0, sigma_u_05, sigma_u_10
    )
    cn = finalize_cn(cn_raw, shape)
    print(f"  cone-native done in {time.time()-t0:.1f}s", flush=True)

    # ----- Family 2: angle / spreadsheet -----
    print("\n=== Family 2: ANGLE (spreadsheet-style) ===", flush=True)
    ang, ang_noise = build_angle_family(
        n_model_grid, pmeta, alpha, d0, sigma_u_05, sigma_u_10
    )

    # ----- Legacy sanity -----
    print("\n=== Legacy sanity ===", flush=True)
    sky = physics_sky_maps(BASE["x"], BASE["f_u"], BASE["alpha"], BASE["eps"])
    geom = build_cone_geometry(L_c, B_c, sky["n_hat"], sky["weights"])
    Z_legacy = cone_average(sky["A_invsq_sky"], geom)
    Z_ref = np.asarray(ref["cone_phys_invsq"], dtype=float)
    d_leg_cn = float(np.nanmax(np.abs(Z_legacy - cn["baseline_invsq"])))
    d_leg_ref = float(np.nanmax(np.abs(Z_legacy - Z_ref)))
    d_cn_ref = float(np.nanmax(np.abs(cn["baseline_invsq"] - Z_ref)))
    d_leg_ang = float(np.nanmax(np.abs(Z_legacy - ang["baseline_invsq"])))
    print(f"  max|legacy − conenative baseline| = {d_leg_cn:.5f}", flush=True)
    print(f"  max|legacy − cone75_ref|          = {d_leg_ref:.5f}", flush=True)
    print(f"  max|conenative baseline − ref|    = {d_cn_ref:.5f}", flush=True)
    print(f"  max|legacy − angle baseline|      = {d_leg_ang:.5f}", flush=True)

    # ----- Cards: cone-native (primary names) -----
    print("\nplotting cards…", flush=True)
    cn_cards = [
        (
            "CARD_cone75_native_baseline_invsq.png",
            cn["baseline_invsq"],
            True,
            rf"θ=75° CONE-NATIVE $A/A_i$ (invsq) — baseline, NO noise"
            "\n"
            rf"$N_{{\mathrm{{samp}}}}={N_SAMP}$ objects/cone | ST+Gauss "
            rf"D/L={BASE['x']}, $f_u$={BASE['f_u']}, α={BASE['alpha']}, d/D={BASE['eps']} "
            rf"| [{VMIN:.2f},{VMAX:.2f}]",
        ),
        (
            "CARD_cone75_native_vel_radial_s05_invsq.png",
            cn["vel_s05_invsq"],
            False,
            rf"θ=75° CONE-NATIVE — radial velocity noise only"
            "\n"
            rf"$\sigma_u=0.05\,H\,d_0$ on each of $N={N_SAMP}$ objects | fixed $d=d_0$ | "
            rf"[{VMIN:.2f},{VMAX:.2f}]",
        ),
        (
            "CARD_cone75_native_vel_radial_s10_invsq.png",
            cn["vel_s10_invsq"],
            False,
            rf"θ=75° CONE-NATIVE — radial velocity noise (louder)"
            "\n"
            rf"$\sigma_u=0.10\,H\,d_0$ | $N={N_SAMP}$ | [{VMIN:.2f},{VMAX:.2f}]",
        ),
        (
            "CARD_cone75_native_dist_s15_invsq.png",
            cn["dist_s15_invsq"],
            False,
            rf"θ=75° CONE-NATIVE — distance scatter only"
            "\n"
            rf"$\sigma_{{d,\mathrm{{rel}}}}=0.15$ per object | $N={N_SAMP}$ | "
            rf"[{VMIN:.2f},{VMAX:.2f}]",
        ),
        (
            "CARD_cone75_native_dist_s25_invsq.png",
            cn["dist_s25_invsq"],
            False,
            rf"θ=75° CONE-NATIVE — distance scatter (louder)"
            "\n"
            rf"$\sigma_{{d,\mathrm{{rel}}}}=0.25$ | $N={N_SAMP}$ | [{VMIN:.2f},{VMAX:.2f}]",
        ),
        (
            "CARD_cone75_native_vel_dist_combined_invsq.png",
            cn["both_s05_s15_invsq"],
            False,
            rf"θ=75° CONE-NATIVE — BOTH velocity + distance noise"
            "\n"
            rf"$\sigma_u=0.05\,H\,d_0$ + $\sigma_{{d,\mathrm{{rel}}}}=0.15$ per object | "
            rf"$N={N_SAMP}$ | [{VMIN:.2f},{VMAX:.2f}]",
        ),
    ]
    for fname, Z, marks, title in cn_cards:
        plot_card(L_c, B_c, Z, title, OUT / fname, show_marks=marks)
        print(f"  → {fname}", flush=True)

    # Angle companion cards
    ang_cards = [
        (
            "CARD_cone75_angle_baseline_invsq.png",
            ang["baseline_invsq"],
            True,
            rf"θ=75° grid ANGLE (spreadsheet-style) — baseline, NO noise"
            "\n"
            rf"one eval per cone-center $\hat{{n}}$ | ST+Gauss D/L={BASE['x']} | "
            rf"[{VMIN:.2f},{VMAX:.2f}]",
        ),
        (
            "CARD_cone75_angle_vel_radial_s05_invsq.png",
            ang["vel_s05_invsq"],
            False,
            rf"θ=75° ANGLE — radial velocity noise only"
            "\n"
            rf"$\sigma_u=0.05\,H\,d_0$ along $\hat{{r}}$ | one eval/angle | "
            rf"[{VMIN:.2f},{VMAX:.2f}]",
        ),
        (
            "CARD_cone75_angle_dist_s15_invsq.png",
            ang["dist_s15_invsq"],
            False,
            rf"θ=75° ANGLE — distance scatter only"
            "\n"
            rf"$\sigma_{{d,\mathrm{{rel}}}}=0.15$ | one eval/angle | [{VMIN:.2f},{VMAX:.2f}]",
        ),
        (
            "CARD_cone75_angle_vel_dist_combined_invsq.png",
            ang["both_s05_s15_invsq"],
            False,
            rf"θ=75° ANGLE — BOTH velocity + distance noise"
            "\n"
            rf"$\sigma_u=0.05\,H\,d_0$ + $\sigma_{{d,\mathrm{{rel}}}}=0.15$ | "
            rf"one eval/angle | [{VMIN:.2f},{VMAX:.2f}]",
        ),
    ]
    for fname, Z, marks, title in ang_cards:
        plot_card(L_c, B_c, Z, title, OUT / fname, show_marks=marks)
        print(f"  → {fname}", flush=True)

    # ----- Strips -----
    def four(maps, prefix_titles=True):
        return [
            {"L": L_c, "B": B_c, "Z": maps["baseline_invsq"],
             "short_title": "baseline\n(no noise)"},
            {"L": L_c, "B": B_c, "Z": maps["vel_s05_invsq"],
             "short_title": r"vel-only" "\n" r"$\sigma_u=0.05 H d_0$"},
            {"L": L_c, "B": B_c, "Z": maps["dist_s15_invsq"],
             "short_title": r"dist-only" "\n" r"$\sigma_d=0.15$"},
            {"L": L_c, "B": B_c, "Z": maps["both_s05_s15_invsq"],
             "short_title": r"vel+dist" "\n" r"combined"},
        ]

    strip_cn = OUT / "STRIP_conenative_baseline_vel_dist_both_invsq.png"
    plot_strip(
        four(cn),
        strip_cn,
        title=(
            rf"CONE-NATIVE θ=75° $A/A_i$ (invsq) — baseline | vel | dist | BOTH"
            rf"  |  $N_{{\mathrm{{samp}}}}={N_SAMP}$  |  magma [{VMIN:.2f},{VMAX:.2f}]"
        ),
    )
    strip_ang = OUT / "STRIP_angle_baseline_vel_dist_both_invsq.png"
    plot_strip(
        four(ang),
        strip_ang,
        title=(
            rf"ANGLE (spreadsheet-style) $A/A_i$ (invsq) — baseline | vel | dist | BOTH"
            rf"  |  magma [{VMIN:.2f},{VMAX:.2f}]"
        ),
    )
    # Back-compat alias for earlier strip name
    strip_alias = OUT / "STRIP_baseline_vel_dist_both_invsq.png"
    plot_strip(
        four(cn),
        strip_alias,
        title=(
            rf"CONE-NATIVE θ=75° $A/A_i$ (invsq) — baseline | vel | dist | BOTH"
            rf"  |  $N={N_SAMP}$  |  magma [{VMIN:.2f},{VMAX:.2f}]"
        ),
    )

    grid_path = OUT / "GRID_conenative_vs_angle_invsq.png"
    plot_compare_grid(
        [
            {"label": "cone-native\n(N=64 in θ=75°)", "maps": four(cn)},
            {"label": "angle\n(1 eval / center)", "maps": four(ang)},
        ],
        grid_path,
        title=(
            rf"Compare: cone-native vs angle  |  baseline | vel | dist | BOTH"
            rf"  |  magma [{VMIN:.2f},{VMAX:.2f}]  |  noise inside $H_{{\mathrm{{obs}}}}$"
        ),
    )
    print(f"strip_cn → {strip_cn.name}", flush=True)
    print(f"strip_ang → {strip_ang.name}", flush=True)
    print(f"grid → {grid_path.name}", flush=True)

    # ----- npz -----
    npz_path = OUT / "cone_param_noise_maps.npz"
    save = dict(
        L_centers=L_c,
        B_centers=B_c,
        theta_deg=THETA_DEG,
        n_samp=N_SAMP,
        x=BASE["x"],
        f_u=BASE["f_u"],
        alpha=BASE["alpha"],
        eps=BASE["eps"],
        d0=d0,
        H=H,
        G=G,
        D=D,
        kappa=pmeta["kappa"],
        migkas_l=MIGKAS_L,
        migkas_b=MIGKAS_B,
        sigma_u_frac_primary=SIGMA_U_FRAC_PRIMARY,
        sigma_u_frac_loud=SIGMA_U_FRAC_LOUD,
        sigma_d_rel_primary=SIGMA_D_REL_PRIMARY,
        sigma_d_rel_loud=SIGMA_D_REL_LOUD,
        sigma_u_primary=sigma_u_05,
        seed_sightlines=SEED_SIGHTLINES,
        seed_cn_vel=SEED_CN_VEL,
        seed_cn_dist=SEED_CN_DIST,
        seed_cn_both_v=SEED_CN_BOTH_V,
        seed_cn_both_d=SEED_CN_BOTH_D,
        seed_ang_vel=SEED_ANG_VEL,
        seed_ang_dist=SEED_ANG_DIST,
        seed_ang_both_v=SEED_ANG_BOTH_V,
        seed_ang_both_d=SEED_ANG_BOTH_D,
        A_legacy_cone75_invsq=Z_legacy,
        A_cone75_ref_invsq=Z_ref,
        legacy_vs_conenative_baseline_max_abs=d_leg_cn,
        legacy_vs_ref_max_abs=d_leg_ref,
        conenative_vs_ref_max_abs=d_cn_ref,
        legacy_vs_angle_baseline_max_abs=d_leg_ang,
    )
    for k, v in cn.items():
        save[f"cn_{k}"] = v
    for k, v in ang.items():
        save[f"ang_{k}"] = v
    # store compact noise (angle full; cn shapes are large — store summaries + seeds only for cn)
    for k, v in ang_noise.items():
        save[f"ang_noise_{k}"] = v
    # cn noise is (2520, 64) — save for reproducibility
    for k, v in cn_noise.items():
        if k == "n_samp":
            continue
        save[f"cn_noise_{k}"] = v

    np.savez_compressed(npz_path, **save)
    print(f"npz → {npz_path}", flush=True)

    # ----- NOTE -----
    def spans_block(maps, tag):
        lines = []
        for k in ["baseline", "vel_s05", "vel_s10", "dist_s15", "dist_s25", "both_s05_s15"]:
            Z = maps[f"{k}_invsq"]
            lines.append(
                f"- **{tag}/{k}**: [{float(np.nanmin(Z)):.4f}, {float(np.nanmax(Z)):.4f}]"
            )
        return "\n".join(lines)

    note = f"""# Cone-param noise: A/A_i with variance INSIDE H_obs

**Created:** {now_local()}

Two families — both kept. Noise is **never** added by shuffling finished map pixels.

## Why ≠ pixel patchwork

| Wrong (patchwork) | This run |
|---|---|
| Build finished $A$, then mix/jitter pixels | Recompute $H_{{\\mathrm{{obs}}}}$ with noisy $u$ or $d$ |
| Noise unrelated to CoE geometry | Velocity noise **along $\\hat{{r}}$ from CoE** |
| Average already-noised sky pixels | Objects (or angles) evaluated natively |

## Family 1 — CONE-NATIVE (primary)

For each cone center $\\hat{{c}}$ on the 72×35 grid from `cone75_maps.npz`:

1. Draw $N_{{\\mathrm{{samp}}}}={N_SAMP}$ sightlines **uniformly in solid angle** within θ={THETA_DEG}° of $\\hat{{c}}$.
2. For each object $\\hat{{n}}_i$: compute model $\\mathbf{{u}}$, optionally add radial velocity noise and/or distance noise, form $H_{{\\mathrm{{obs}}}}$ and $A_{{\\mathrm{{raw}}}}$.
3. Cone $A$ = mean of $A_{{\\mathrm{{raw}}}}$ over the $N_{{\\mathrm{{samp}}}}$ objects.
4. Mean-renormalize the cone map → $A/A_i$.

Shared sightline draws: seed `{SEED_SIGHTLINES}`. Noise seeds: vel={SEED_CN_VEL},
dist={SEED_CN_DIST}, both-v={SEED_CN_BOTH_V}, both-d={SEED_CN_BOTH_D}.

Cards: `CARD_cone75_native_*.png`
Strip: `STRIP_conenative_baseline_vel_dist_both_invsq.png`

## Family 2 — ANGLE (spreadsheet-style companion)

One evaluation per cone-center direction $\\hat{{n}}$ (the spreadsheet idea:
add random to velocity/distance **in the calculation** at that angle). Same
physics dials and σ choices. Independent seeds (100+ offsets vs cone-native).

Cards: `CARD_cone75_angle_*.png`
Strip: `STRIP_angle_baseline_vel_dist_both_invsq.png`

## Variants (each family)

| map | rule |
|-----|------|
| **baseline** | $d=d_0=\\varepsilon D$, pure model $u$ |
| **vel-only** | $d=d_0$; $\\mathbf{{u}}_r \\leftarrow \\mathbf{{u}}_r + \\xi\\sigma_u\\hat{{r}}$ |
| **dist-only** | $d_i=d_0(1+\\eta\\sigma_{{d,\\mathrm{{rel}}}})$, clip $>0.02\\,d_0$; model $u$ |
| **BOTH** | vel + dist noise together on the same objects/angles |

### Rationalized σ

- **$\\sigma_u = 0.05\\,H\\,d_0$** (primary; also 0.10 louder): fraction of shell Hubble
  velocity $V_H=H d_0$. Mild few-percent hit on $H_{{\\mathrm{{obs}}}}$; order-of-magnitude
  akin to few×10²–10³ km/s if mapped to local-cluster scales (illustrative, not
  a cosmographic conversion). Noise is **radial from CoE** ($\\hat{{r}}$).
- **$\\sigma_{{d,\\mathrm{{rel}}}} = 0.15$** (primary; also 0.25): ~10–20% depth/measurement
  mix is plausible for Migkas-like wide-$z$ samples within a large cone.

## Physics dials
- twoscale local-match, ρ(D)=Ω_m ρ_crit
- D/L={BASE['x']}, f_u={BASE['f_u']}, α={BASE['alpha']}, ε={BASE['eps']}
- κ ≈ {pmeta['kappa']:.6f}; units H=G=D=1
- CoE $+x$ → Migkas $(\\ell,b)=({MIGKAS_L}^\\circ,{MIGKAS_B}^\\circ)$

## Legacy sanity vs `cone75_maps.npz`

Legacy = fine full-sky $H_{{\\mathrm{{obs}}}}$→$A$→θ=75° solid-angle cone average.

| comparison | max\\|Δ\\| |
|---|---|
| legacy vs **cone-native** baseline | {d_leg_cn:.5f} |
| legacy vs cone75_ref | {d_leg_ref:.5f} |
| cone-native baseline vs ref | {d_cn_ref:.5f} |
| legacy vs **angle** baseline | {d_leg_ang:.5f} |

Cone-native no-noise should be **close** to legacy/ref (same cone idea; Monte Carlo
vs pixel integral → small residual). Angle baseline is a different operator (point
sample) — larger Δ expected.

## Plotting
- Migkas Galactic: ℓ LEFT (`lon_disp=−ℓ_wrapped`)
- magma [{VMIN:.2f}, {VMAX:.2f}], cards ~{CARD_PX[0]}×{CARD_PX[1]}
- Compare grid: `GRID_conenative_vs_angle_invsq.png`

## A/A_i spans (invsq, mean-renormed)

### Cone-native
{spans_block(cn, "cn")}

### Angle
{spans_block(ang, "ang")}

## Files
- `run_cone_param_noise.py`, `cone_param_noise_maps.npz`, `manifest.json`, this NOTE
"""
    (OUT / "NOTE.md").write_text(note)

    man = {
        "created_local": now_local(),
        "families": ["conenative", "angle"],
        "n_samp": N_SAMP,
        "theta_deg": THETA_DEG,
        "baseline": BASE,
        "kappa": float(pmeta["kappa"]),
        "sigma_u_frac_primary": SIGMA_U_FRAC_PRIMARY,
        "sigma_d_rel_primary": SIGMA_D_REL_PRIMARY,
        "legacy": {
            "vs_conenative_baseline_max_abs": d_leg_cn,
            "vs_ref_max_abs": d_leg_ref,
            "conenative_vs_ref_max_abs": d_cn_ref,
            "vs_angle_baseline_max_abs": d_leg_ang,
        },
        "spans_conenative": {
            k: [float(np.nanmin(cn[f"{k}_invsq"])), float(np.nanmax(cn[f"{k}_invsq"]))]
            for k in ["baseline", "vel_s05", "vel_s10", "dist_s15", "dist_s25", "both_s05_s15"]
        },
        "spans_angle": {
            k: [float(np.nanmin(ang[f"{k}_invsq"])), float(np.nanmax(ang[f"{k}_invsq"]))]
            for k in ["baseline", "vel_s05", "vel_s10", "dist_s15", "dist_s25", "both_s05_s15"]
        },
        "strip_conenative": str(strip_cn),
        "strip_angle": str(strip_ang),
        "grid": str(grid_path),
        "npz": str(npz_path),
        "elapsed_s": round(time.time() - t_all, 2),
    }
    (OUT / "manifest.json").write_text(json.dumps(man, indent=2))
    print(f"\nDone in {time.time()-t_all:.1f}s", flush=True)
    return man


if __name__ == "__main__":
    main()

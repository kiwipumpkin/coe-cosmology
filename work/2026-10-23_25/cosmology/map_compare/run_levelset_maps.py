#!/usr/bin/env python3
"""
NO-NOISE θ=75° cone A/A_i Mollweide cards across level-set dials.

Family A: baseline ST+Gauss physics fixed; Amp/φ diagnostic only.
Family B: contour level-set peaks (profile dials change).

Plotting: Migkas Galactic (ℓ increases LEFT: lon_disp = −ℓ_wrapped),
magma [0.85, 1.20], 1400×900 cards, celestial equator if astropy available.
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
from mpl_toolkits.axes_grid1 import make_axes_locatable

ROOT = Path("/workspace/cosmology")
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "density_local_match"))
sys.path.insert(0, str(ROOT / "density_local_match" / "r2_maximize"))

from run_density_local_match import (  # noqa: E402
    params_twoscale_local_match,
    peculiar_u,
    m_enc_at_r,
    D,
    H,
)
from forward_model import make_latlon_grid  # noqa: E402

OUT = Path(__file__).resolve().parent / "levelset_maps"
FULLSKY = ROOT / "fullsky_AB" / "fullsky_AB_maps.npz"
CONE_REF = ROOT / "map_compare" / "cone75_maps.npz"

THETA_DEG = 75.0
VMIN, VMAX = 0.85, 1.20
OMEGA_L = 2.0 * np.pi / 360.0
MIGKAS_L, MIGKAS_B = 280.0, -20.0
CARD_PX = (1400, 900)
DPI = 100  # 14"×9" @ 100 dpi = 1400×900

BASE = dict(x=5.9, f_u=0.94, alpha=1.0, eps=0.1)

try:
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    HAS_ASTROPY = True
except Exception:
    HAS_ASTROPY = False


def now_local() -> str:
    return datetime.now(ZoneInfo("America/Chicago")).strftime(
        "%Y-%m-%d %I:%M %p %Z (UTC%z)"
    )


def wrap_lon_deg(lon: np.ndarray) -> np.ndarray:
    """Wrap to (−180, 180]."""
    return ((np.asarray(lon, dtype=float) + 180.0) % 360.0) - 180.0


def lon_disp_from_ell(ell_deg: np.ndarray) -> np.ndarray:
    """Migkas Galactic: ℓ increases LEFT → plot lon_disp = −ℓ_wrapped."""
    return -wrap_lon_deg(ell_deg)


def rotation_matrix_x_to_sky(l_deg: float, b_deg: float) -> np.ndarray:
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)
    n_gal = np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])
    ez = np.array([0.0, 0.0, 1.0])
    ytmp = np.cross(ez, n_gal)
    if np.linalg.norm(ytmp) < 1e-8:
        ytmp = np.cross(np.array([0.0, 1.0, 0.0]), n_gal)
    ey = ytmp / np.linalg.norm(ytmp)
    ez_new = np.cross(n_gal, ey)
    ez_new /= np.linalg.norm(ez_new)
    return np.column_stack([n_gal, ey, ez_new])


def rotate_map_to_sky(values, Lon, Lat, R, n_lon, n_lat):
    nx = np.cos(Lat) * np.cos(Lon)
    ny = np.cos(Lat) * np.sin(Lon)
    nz = np.sin(Lat)
    n_sky = np.stack([nx.ravel(), ny.ravel(), nz.ravel()], axis=-1)
    n_model = (R.T @ n_sky.T).T
    lon_m = np.arctan2(n_model[:, 1], n_model[:, 0]) % (2.0 * np.pi)
    lat_m = np.arcsin(np.clip(n_model[:, 2], -1.0, 1.0))
    V = values.reshape(n_lat, n_lon)
    lon_1d = np.linspace(0.0, 2.0 * np.pi, n_lon, endpoint=False)
    lat_1d = np.linspace(-0.5 * np.pi + 1e-3, 0.5 * np.pi - 1e-3, n_lat)
    fi = lon_m / (2.0 * np.pi) * n_lon
    lat0, lat1 = lat_1d[0], lat_1d[-1]
    fj = (lat_m - lat0) / (lat1 - lat0) * (n_lat - 1)
    i0 = np.floor(fi).astype(int) % n_lon
    i1 = (i0 + 1) % n_lon
    j0 = np.clip(np.floor(fj).astype(int), 0, n_lat - 2)
    j1 = j0 + 1
    wi = fi - np.floor(fi)
    wj = np.clip(fj - np.floor(fj), 0.0, 1.0)
    out = (
        (1 - wi) * (1 - wj) * V[j0, i0]
        + wi * (1 - wj) * V[j0, i1]
        + (1 - wi) * wj * V[j1, i0]
        + wi * wj * V[j1, i1]
    )
    return out.reshape(n_lat, n_lon)


def apply_sine_A(AF, l_deg, amp_a, phi):
    """A' = (1 − Amp_A · sin(2πℓ/360+φ)) · A_physics"""
    s = np.sin(OMEGA_L * np.asarray(l_deg, dtype=float) + phi)
    return (1.0 - amp_a * s) * AF


def H_obs_sphere(pmeta, eps, alpha, n_hat):
    d = eps * D
    D_vec = np.array([D, 0.0, 0.0])
    r_vec = D_vec + d * n_hat
    u_r = peculiar_u(r_vec, pmeta, alpha=alpha)
    u_D = peculiar_u(D_vec.reshape(1, 3), pmeta, alpha=alpha)[0]
    delta_v = np.sum((u_r - u_D) * n_hat, axis=-1)
    H_obs = H + delta_v / d
    m_src = m_enc_at_r(np.linalg.norm(r_vec, axis=-1), pmeta)
    return H_obs, m_src


def physics_sky_maps(x, f_u, alpha, eps, n_lon=360, n_lat=181):
    """Return A_invsq_sky, A_hratio_sky, l_gal_deg, weights, Lon, Lat (radians)."""
    pmeta = params_twoscale_local_match(x, f_u)
    n_hat, lon_flat, lat_flat, w = make_latlon_grid(n_lon, n_lat)
    Lon = lon_flat.reshape(n_lat, n_lon)
    Lat = lat_flat.reshape(n_lat, n_lon)
    H_obs, m_src = H_obs_sphere(pmeta, eps, alpha, n_hat)
    H_safe = np.maximum(H_obs, 1e-6 * H)
    A_hr_raw = H_safe / H
    A_inv_raw = 1.0 / (H_safe / H) ** 2

    def mean_renorm(xarr):
        m = float(np.average(xarr, weights=w))
        return (xarr / m).reshape(n_lat, n_lon)

    A_hr = mean_renorm(A_hr_raw)
    A_inv = mean_renorm(A_inv_raw)
    R = rotation_matrix_x_to_sky(MIGKAS_L, MIGKAS_B)
    A_hr_sky = rotate_map_to_sky(A_hr, Lon, Lat, R, n_lon, n_lat)
    A_inv_sky = rotate_map_to_sky(A_inv, Lon, Lat, R, n_lon, n_lat)
    l_gal = (np.rad2deg(Lon) % 360.0)  # sky Lon ≡ Galactic ℓ
    w2 = w.reshape(n_lat, n_lon)
    return {
        "A_invsq_sky": A_inv_sky,
        "A_hratio_sky": A_hr_sky,
        "l_gal_deg": l_gal,
        "weights": w2,
        "Lon": Lon,
        "Lat": Lat,
        "pmeta": pmeta,
        "n_hat": n_hat.reshape(n_lat, n_lon, 3),
    }


def build_cone_geometry(L_centers, B_centers, n_hat_pix, weights):
    """Precompute unit vectors for centers and flatten pixel geometry."""
    Lc, Bc = np.meshgrid(L_centers, B_centers)  # (nB, nL)
    l = np.deg2rad(Lc.ravel())
    b = np.deg2rad(Bc.ravel())
    c = np.stack(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], axis=-1
    )
    n_flat = n_hat_pix.reshape(-1, 3)
    w_flat = weights.ravel()
    cos_th = np.cos(np.deg2rad(THETA_DEG))
    return c, n_flat, w_flat, cos_th, Lc.shape


def cone_average(values, geom):
    c, n_flat, w_flat, cos_th, shape = geom
    v_flat = values.ravel()
    # chunk centers to limit memory
    out = np.empty(c.shape[0], dtype=float)
    chunk = 200
    for i0 in range(0, c.shape[0], chunk):
        i1 = min(i0 + chunk, c.shape[0])
        dots = c[i0:i1] @ n_flat.T  # (nc, npix)
        for j in range(i1 - i0):
            m = dots[j] >= cos_th
            ww = w_flat[m]
            sw = ww.sum()
            out[i0 + j] = np.sum(ww * v_flat[m]) / sw if sw > 0 else np.nan
    Z = out.reshape(shape)
    return Z / np.nanmean(Z)


def celestial_equator_lb(n=720):
    if not HAS_ASTROPY:
        return None, None
    ra = np.linspace(0.0, 360.0, n, endpoint=False)
    sc = SkyCoord(ra=ra * u.deg, dec=np.zeros_like(ra) * u.deg, frame="icrs")
    gal = sc.galactic
    return np.asarray(gal.l.deg), np.asarray(gal.b.deg)


def style():
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "font.size": 11,
            "axes.titlesize": 13,
            "mathtext.fontset": "dejavusans",
        }
    )


def plot_card(L_c, B_c, Z, title, path, show_marks=False):
    """1400×900 Migkas-convention Mollweide card, magma [0.85,1.20]."""
    style()
    # cell edges for pcolormesh
    dL = float(np.median(np.diff(L_c))) if len(L_c) > 1 else 5.0
    dB = float(np.median(np.diff(B_c))) if len(B_c) > 1 else 5.0
    L_edge = np.concatenate([L_c - dL / 2, [L_c[-1] + dL / 2]])
    B_edge = np.concatenate([B_c - dB / 2, [B_c[-1] + dB / 2]])
    # wrap L to continuous for display via lon_disp
    # Build edge mesh in display coords
    # For each column of Z corresponding to L_c[j], lon_disp = -wrap(L)
    Lon_c, Lat_c = np.meshgrid(L_c, B_c)
    lon_disp = np.deg2rad(lon_disp_from_ell(Lon_c))
    lat = np.deg2rad(Lat_c)

    # Sort columns by lon_disp so pcolormesh is happy
    order = np.argsort(lon_disp[0])
    lon_s = lon_disp[:, order]
    lat_s = lat[:, order]
    Z_s = Z[:, order]

    # Build edges from sorted centers
    lon_centers = lon_s[0]
    # midpoint edges
    mid = 0.5 * (lon_centers[:-1] + lon_centers[1:])
    # extrapolate ends
    lon_edges_1d = np.concatenate(
        [
            [lon_centers[0] - (mid[0] - lon_centers[0])],
            mid,
            [lon_centers[-1] + (lon_centers[-1] - mid[-1])],
        ]
    )
    lat_edges_1d = np.deg2rad(B_edge)
    Lon_e, Lat_e = np.meshgrid(lon_edges_1d, lat_edges_1d)

    fig_w, fig_h = CARD_PX[0] / DPI, CARD_PX[1] / DPI
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=DPI)
    ax = fig.add_subplot(111, projection="mollweide")
    pcm = ax.pcolormesh(
        Lon_e,
        Lat_e,
        Z_s,
        shading="flat",
        cmap="magma",
        norm=Normalize(VMIN, VMAX),
        rasterized=True,
    )
    ax.grid(True, linestyle=":", alpha=0.45, color="white", linewidth=0.6)

    # Celestial equator
    ell_eq, b_eq = celestial_equator_lb()
    if ell_eq is not None:
        ld = np.deg2rad(lon_disp_from_ell(ell_eq))
        bd = np.deg2rad(b_eq)
        # break at wrap jumps
        jumps = np.where(np.abs(np.diff(ld)) > np.pi / 2)[0]
        starts = np.concatenate([[0], jumps + 1])
        ends = np.concatenate([jumps + 1, [len(ld)]])
        for s, e in zip(starts, ends):
            if e - s < 2:
                continue
            ax.plot(ld[s:e], bd[s:e], "k-", lw=1.2, alpha=0.85, zorder=5)

    # Longitude tick labels: 90 (left), 0 (center), 270 (right)
    # In display coords: lon_disp(90)=-90°=left, lon_disp(0)=0, lon_disp(270)=+90°=right
    ax.set_xticks(np.deg2rad([-90, 0, 90]))
    ax.set_xticklabels([r"$90^\circ$", r"$0^\circ$", r"$270^\circ$"], fontsize=11)
    ax.set_yticks(np.deg2rad([-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75]))
    ax.tick_params(axis="y", labelsize=8)

    if show_marks:
        # Migkas sample min approx
        for ell, bb, marker, color, lab in [
            (280.0, -20.0, "*", "white", r"Migkas $(280^\circ,-20^\circ)$"),
        ]:
            ax.plot(
                np.deg2rad(lon_disp_from_ell(ell)),
                np.deg2rad(bb),
                marker=marker,
                color=color,
                ms=10,
                mew=0.8,
                mec="k",
                zorder=6,
                label=lab,
            )
        ax.legend(loc="lower left", fontsize=8, framealpha=0.7)

    ax.set_title(title, pad=14, fontsize=12)

    # Colorbar with numbered ticks
    cbar = fig.colorbar(
        pcm, ax=ax, orientation="vertical", pad=0.02, fraction=0.035, aspect=28
    )
    ticks = np.linspace(VMIN, VMAX, 8)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t:.2f}" for t in ticks])
    cbar.set_label(r"$A/A_i$", fontsize=11)
    cbar.ax.tick_params(labelsize=9)

    # Generous margins so full Mollweide visible
    fig.subplots_adjust(left=0.04, right=0.92, bottom=0.06, top=0.90)
    fig.savefig(path, dpi=DPI, facecolor="white")
    plt.close(fig)
    # enforce exact pixel size
    try:
        from PIL import Image

        im = Image.open(path)
        if im.size != CARD_PX:
            im = im.resize(CARD_PX, Image.Resampling.LANCZOS)
            im.save(path)
    except Exception:
        pass


def plot_grid(entries, path, ncols=3, title="Amp×sin family — θ=75° cone A/A_i"):
    """Summary grid with shared colorbar; equal subplot sizes."""
    style()
    n = len(entries)
    nrows = int(np.ceil(n / ncols))
    fig_w = 4.6 * ncols + 1.2
    fig_h = 2.7 * nrows + 0.9
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(fig_w, fig_h),
        subplot_kw={"projection": "mollweide"},
        constrained_layout=False,
    )
    axes = np.atleast_1d(axes).ravel()
    pcm = None
    for ax, ent in zip(axes, entries):
        L_c, B_c, Z = ent["L"], ent["B"], ent["Z"]
        Lon_c, Lat_c = np.meshgrid(L_c, B_c)
        lon_disp = np.deg2rad(lon_disp_from_ell(Lon_c))
        lat = np.deg2rad(Lat_c)
        order = np.argsort(lon_disp[0])
        lon_s = lon_disp[:, order]
        Z_s = Z[:, order]
        lat_s = lat[:, order]
        lon_centers = lon_s[0]
        mid = 0.5 * (lon_centers[:-1] + lon_centers[1:])
        lon_edges = np.concatenate(
            [
                [lon_centers[0] - (mid[0] - lon_centers[0])],
                mid,
                [lon_centers[-1] + (lon_centers[-1] - mid[-1])],
            ]
        )
        dB = float(np.median(np.diff(B_c)))
        B_edge = np.concatenate([B_c - dB / 2, [B_c[-1] + dB / 2]])
        Lon_e, Lat_e = np.meshgrid(lon_edges, np.deg2rad(B_edge))
        pcm = ax.pcolormesh(
            Lon_e,
            Lat_e,
            Z_s,
            shading="flat",
            cmap="magma",
            norm=Normalize(VMIN, VMAX),
            rasterized=True,
        )
        ax.grid(True, linestyle=":", alpha=0.35, color="white", lw=0.4)
        ax.set_xticks(np.deg2rad([-90, 0, 90]))
        ax.set_xticklabels(["90", "0", "270"], fontsize=7)
        ax.set_yticks([])
        ax.set_title(ent["short_title"], fontsize=8, pad=8)
    for ax in axes[n:]:
        ax.set_visible(False)

    fig.subplots_adjust(
        left=0.03, right=0.90, bottom=0.04, top=0.90, wspace=0.08, hspace=0.28
    )
    fig.suptitle(title + f"  |  magma [{VMIN:.2f}, {VMAX:.2f}]  |  NO noise", fontsize=12)
    if pcm is not None:
        cax = fig.add_axes([0.92, 0.12, 0.015, 0.7])
        cbar = fig.colorbar(pcm, cax=cax)
        cbar.set_label(r"$A/A_i$")
        ticks = np.linspace(VMIN, VMAX, 8)
        cbar.set_ticks(ticks)
    fig.savefig(path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def load_baseline_sky():
    fs = np.load(FULLSKY)
    keys = list(fs.keys())
    inv = fs["A_invsq_sky"]
    hr = fs["A_hratio_sky"]
    if "l_gal_sky_deg" in keys:
        l_gal = fs["l_gal_sky_deg"]
    elif "l_gal_deg" in keys:
        l_gal = fs["l_gal_deg"]
    else:
        Lon = fs["lon"]
        l_gal = (np.rad2deg(Lon) % 360.0) if np.nanmax(np.abs(Lon)) <= 7 else (Lon % 360.0)
    w = fs["weights"] if "weights" in keys else fs["weights"]
    Lon = fs["lon"]
    Lat = fs["lat"]
    if "n_hat" in keys:
        n_hat = fs["n_hat"]
    else:
        n_hat = np.stack(
            [
                np.cos(Lat) * np.cos(Lon),
                np.cos(Lat) * np.sin(Lon),
                np.sin(Lat),
            ],
            axis=-1,
        )
    return {
        "A_invsq_sky": inv,
        "A_hratio_sky": hr,
        "l_gal_deg": np.asarray(l_gal, dtype=float),
        "weights": w,
        "Lon": Lon,
        "Lat": Lat,
        "n_hat": n_hat,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t_all = time.time()
    print(f"[{now_local()}] starting level-set cone maps → {OUT}", flush=True)

    ref = np.load(CONE_REF)
    L_c = np.asarray(ref["L_centers"], dtype=float)
    B_c = np.asarray(ref["B_centers"], dtype=float)
    print(f"grid: L={L_c.shape} B={B_c.shape} θ={THETA_DEG}°", flush=True)

    sky0 = load_baseline_sky()
    geom = build_cone_geometry(L_c, B_c, sky0["n_hat"], sky0["weights"])
    print("cone geometry ready", flush=True)

    # Verify phys match
    Z_check = cone_average(sky0["A_invsq_sky"], geom)
    print(
        f"baseline phys invsq cone range [{Z_check.min():.4f}, {Z_check.max():.4f}] "
        f"(ref [{ref['cone_phys_invsq'].min():.4f}, {ref['cone_phys_invsq'].max():.4f}])",
        flush=True,
    )

    # ---------- Family A: Amp/φ variants (baseline physics) ----------
    phi_sheet = 0.35
    phi_local = float(np.deg2rad(212.81767955801106))
    amp_specs = [
        # (nn, shortname, Amp, phi_rad, title extras)
        ("01", "amp0_physics", 0.0, 0.0, "Amp=0 (physics only)"),
        ("02", "amp_p003_phi035", +0.03, phi_sheet, r"Amp=+0.03, φ=0.35 rad"),
        ("03", "amp_m003_phi035", -0.03, phi_sheet, r"Amp=−0.03, φ=0.35 rad"),
        ("04", "amp_p006_phi035", +0.06, phi_sheet, r"Amp=+0.06, φ=0.35 rad"),
        ("05", "amp_m006_phi035_sheet", -0.06, phi_sheet, r"Amp=−0.06, φ=0.35 (sheet)"),
        ("06", "amp_p010_phi035", +0.10, phi_sheet, r"Amp=+0.10, φ=0.35 rad"),
        ("07", "amp_m010_phi035", -0.10, phi_sheet, r"Amp=−0.10, φ=0.35 rad"),
        (
            "08",
            "amp_p007_phi213_local",
            +0.07,
            phi_local,
            r"Amp=+0.07, φ≈212.82° (local-density)",
        ),
    ]

    manifest = {
        "created_local": now_local(),
        "theta_deg": THETA_DEG,
        "af_display": "invsq",
        "colorbar": [VMIN, VMAX],
        "cmap": "magma",
        "noise": False,
        "plot_convention": "Migkas Galactic: lon_disp=-ell_wrapped (ell increases LEFT)",
        "baseline_physics": dict(BASE),
        "maps": [],
        "failures": [],
    }
    npz_data = {"L_centers": L_c, "B_centers": B_c, "theta_deg": THETA_DEG}
    amp_grid_entries = []

    for nn, short, amp, phi, label in amp_specs:
        t0 = time.time()
        A_pix = apply_sine_A(sky0["A_invsq_sky"], sky0["l_gal_deg"], amp, phi)
        A_hr_pix = apply_sine_A(sky0["A_hratio_sky"], sky0["l_gal_deg"], amp, phi)
        Z = cone_average(A_pix, geom)
        Z_hr = cone_average(A_hr_pix, geom)
        title = (
            rf"θ=75° cone $A/A_i$ (invsq) — {label}"
            "\n"
            rf"baseline physics D/L=5.9, $f_u$=0.94, α=1, d/D=0.1 | NO noise | [{VMIN:.2f},{VMAX:.2f}]"
        )
        card = OUT / f"map_{nn}_{short}.png"
        plot_card(L_c, B_c, Z, title, card, show_marks=(amp == 0.0))
        key = f"A_{short}"
        npz_data[key] = Z
        npz_data[f"{key}_hratio"] = Z_hr
        entry = {
            "nn": nn,
            "name": short,
            "family": "A_amp_phi",
            "Amp": amp,
            "phi_rad": phi,
            "phi_deg": float(np.rad2deg(phi)) if amp != 0 else 0.0,
            "D_over_L": BASE["x"],
            "f_u": BASE["f_u"],
            "alpha": BASE["alpha"],
            "d_over_D": BASE["eps"],
            "af": "invsq",
            "A_min": float(np.nanmin(Z)),
            "A_max": float(np.nanmax(Z)),
            "A_mean": float(np.nanmean(Z)),
            "hratio_min": float(np.nanmin(Z_hr)),
            "hratio_max": float(np.nanmax(Z_hr)),
            "card": str(card),
            "kind": "physics" if amp == 0 else "Amp_sin_diagnostic",
            "dt_s": round(time.time() - t0, 2),
        }
        manifest["maps"].append(entry)
        amp_grid_entries.append(
            {
                "L": L_c,
                "B": B_c,
                "Z": Z,
                "short_title": f"{nn} {label}",
            }
        )
        print(
            f"  [{nn}] {short}: Amp={amp:+.3f} φ={entry['phi_deg']:.1f}° "
            f"A∈[{entry['A_min']:.4f},{entry['A_max']:.4f}] → {card.name} ({entry['dt_s']}s)",
            flush=True,
        )

    grid_amp = OUT / "GRID_amp_family.png"
    plot_grid(
        amp_grid_entries,
        grid_amp,
        ncols=3,
        title="A) Amp×sin diagnostic family (baseline ST+Gauss fixed)",
    )
    print(f"Amp grid → {grid_amp}", flush=True)

    # ---------- Family B: profile peaks ----------
    b_specs = [
        (
            "09",
            "peak_baseline_rec",
            dict(x=5.9, f_u=0.94, alpha=1.0, eps=0.1),
            +0.07,
            phi_local,
            "Baseline recommended + Amp=+0.07, φ≈212.82°",
        ),
        (
            "10",
            "peak_DL_fu",
            dict(x=5.69, f_u=0.98, alpha=1.0, eps=0.1),
            +0.02,
            float(np.deg2rad(217.0)),
            "D/L–f_u peak: D/L=5.69, f_u=0.98 + Amp=+0.02, φ=217°",
        ),
        (
            "11",
            "peak_alpha_eps",
            dict(x=5.9, f_u=0.94, alpha=0.62, eps=0.22),
            0.0,
            0.0,
            "α–d/D peak: α=0.62, d/D=0.22 (physics; no Amp in NOTE)",
        ),
        (
            "12",
            "peak_DL_alpha",
            dict(x=5.69, f_u=0.94, alpha=0.50, eps=0.1),
            0.0,
            0.0,
            "D/L–α peak: D/L=5.69, α=0.50 (physics; no Amp in NOTE)",
        ),
    ]

    peak_grid_entries = []
    for nn, short, prof, amp, phi, label in b_specs:
        t0 = time.time()
        try:
            # Reuse baseline sky if profile matches BASE (save time / exact match)
            if (
                abs(prof["x"] - BASE["x"]) < 1e-9
                and abs(prof["f_u"] - BASE["f_u"]) < 1e-9
                and abs(prof["alpha"] - BASE["alpha"]) < 1e-9
                and abs(prof["eps"] - BASE["eps"]) < 1e-9
            ):
                sky = sky0
                pmeta = params_twoscale_local_match(prof["x"], prof["f_u"])
            else:
                print(f"  recomputing physics for {short} …", flush=True)
                sky = physics_sky_maps(
                    prof["x"], prof["f_u"], prof["alpha"], prof["eps"]
                )
                pmeta = sky["pmeta"]
                # rebuild geom if n_hat differs (same grid actually)
                geom_b = build_cone_geometry(L_c, B_c, sky["n_hat"], sky["weights"])
            geom_use = geom if sky is sky0 else geom_b

            A_pix = apply_sine_A(sky["A_invsq_sky"], sky["l_gal_deg"], amp, phi)
            A_hr_pix = apply_sine_A(sky["A_hratio_sky"], sky["l_gal_deg"], amp, phi)
            Z = cone_average(A_pix, geom_use)
            Z_hr = cone_average(A_hr_pix, geom_use)
            title = (
                rf"θ=75° cone $A/A_i$ (invsq) — {label}"
                "\n"
                rf"D/L={prof['x']}, $f_u$={prof['f_u']}, α={prof['alpha']}, d/D={prof['eps']} "
                rf"| κ≈{pmeta['kappa']:.3f} | NO noise | [{VMIN:.2f},{VMAX:.2f}]"
            )
            card = OUT / f"map_{nn}_{short}.png"
            plot_card(L_c, B_c, Z, title, card)
            key = f"A_{short}"
            npz_data[key] = Z
            npz_data[f"{key}_hratio"] = Z_hr
            entry = {
                "nn": nn,
                "name": short,
                "family": "B_profile_peak",
                "Amp": amp,
                "phi_rad": phi,
                "phi_deg": float(np.rad2deg(phi)) if amp != 0 else 0.0,
                "D_over_L": prof["x"],
                "f_u": prof["f_u"],
                "alpha": prof["alpha"],
                "d_over_D": prof["eps"],
                "kappa": float(pmeta["kappa"]),
                "af": "invsq",
                "A_min": float(np.nanmin(Z)),
                "A_max": float(np.nanmax(Z)),
                "A_mean": float(np.nanmean(Z)),
                "hratio_min": float(np.nanmin(Z_hr)),
                "hratio_max": float(np.nanmax(Z_hr)),
                "card": str(card),
                "kind": "physics" if amp == 0 else "physics+Amp_sin",
                "dt_s": round(time.time() - t0, 2),
            }
            manifest["maps"].append(entry)
            peak_grid_entries.append(
                {"L": L_c, "B": B_c, "Z": Z, "short_title": f"{nn} {short}"}
            )
            print(
                f"  [{nn}] {short}: A∈[{entry['A_min']:.4f},{entry['A_max']:.4f}] "
                f"→ {card.name} ({entry['dt_s']}s)",
                flush=True,
            )
        except Exception as e:
            msg = f"{nn} {short}: {type(e).__name__}: {e}"
            print("FAIL", msg, flush=True)
            manifest["failures"].append(msg)

    if peak_grid_entries:
        grid_peak = OUT / "GRID_profile_peaks.png"
        plot_grid(
            peak_grid_entries,
            grid_peak,
            ncols=2,
            title="B) Contour level-set peaks (profile dials)",
        )
        print(f"Peak grid → {grid_peak}", flush=True)

    npz_path = OUT / "levelset_cone75_maps.npz"
    np.savez_compressed(npz_path, **npz_data)
    man_path = OUT / "manifest.json"
    manifest["npz"] = str(npz_path)
    manifest["grid_amp_family"] = str(grid_amp)
    if peak_grid_entries:
        manifest["grid_profile_peaks"] = str(OUT / "GRID_profile_peaks.png")
    man_path.write_text(json.dumps(manifest, indent=2))

    note = f"""# Level-set θ=75° cone A/A_i maps (NO noise)

**Created:** {manifest['created_local']}

## What is physics vs Amp×sin diagnostic

### Physics (ST+Gauss two-scale)
- Profile: twoscale with ρ(D)=Ω_m ρ_crit=0.3 ρ_crit
- Dials: **D/L**, **f_u**, **α**, **d/D**
- AF for these cards: **invsq** — A_raw ∝ 1/(H_obs/H)², mean-renormed on the sky,
  then θ=75° solid-angle cone average, then mean-renorm → A/A_i
- Also stored: **hratio** companion maps in the npz (A_raw = H_obs/H)

### Amp×sin diagnostic (NOT first-principles)
Applied on the sky as a function of Galactic longitude ℓ:

```
A' = (1 − Amp_A · sin(2π ℓ/360 + φ)) · A_physics
```

then cone-averaged and mean-normalized. This is the same sheet/local-density
diagnostic used in `cone_scan_75` / `fullsky_AB` — **not** derived from Gauss.

## Families

### A) Amp×sin diagnostic (baseline physics fixed)
Baseline: D/L=5.9, f_u=0.94, α=1, d/D=0.1 (recommended profile).
Only Amp and φ change — matches the level-set CURVE family on A vs ℓ.

| # | name | Amp | φ |
|---|------|-----|---|
| 01 | amp0_physics | 0 | — |
| 02–07 | ±0.03, ±0.06, ±0.10 | at φ=0.35 rad (sheet phase) |
| 08 | local-density | +0.07 | ≈212.82° |

### B) Contour level-set peaks (profile dials)
| # | peak | dials | Amp/φ |
|---|------|-------|-------|
| 09 | baseline recommended | D/L=5.9, f_u=0.94 | +0.07 / 212.82° |
| 10 | D/L vs f_u | D/L=5.69, f_u=0.98 | +0.02 / 217° (from NOTE) |
| 11 | α vs d/D | α=0.62, d/D=0.22 | physics only (no Amp in NOTE) |
| 12 | D/L vs α | D/L=5.69, α=0.50 | physics only (no Amp in NOTE) |

## Plotting convention
- Migkas Galactic: **ℓ increases LEFT** (`lon_disp = −ℓ_wrapped`)
- magma colorbar fixed **[0.85, 1.20]**
- Cards **1400×900**, full Mollweide visible, ℓ labels 90 / 0 / 270
- Celestial equator via astropy SkyCoord (if available)
- **NO noise** on any map

## Files
- Individual cards: `map_NN_<shortname>.png`
- `GRID_amp_family.png`, `GRID_profile_peaks.png`
- `levelset_cone75_maps.npz`, `manifest.json`

Preferred axis for sky rotation: (ℓ, b) = (280°, −20°).
"""
    (OUT / "NOTE.md").write_text(note)
    print(f"\nDone in {time.time()-t_all:.1f}s", flush=True)
    print(f"manifest → {man_path}", flush=True)
    print(f"npz → {npz_path}", flush=True)


if __name__ == "__main__":
    main()

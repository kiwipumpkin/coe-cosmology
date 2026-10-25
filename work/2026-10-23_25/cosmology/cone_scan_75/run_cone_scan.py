#!/usr/bin/env python3
"""
θ=75° cone-scan of full-sky A/B maps vs Migkas Fig.6 strips.

PRIMARY: pure ST+Gauss physics averages (A_hratio, A_invsq, B=AN mass-proxy).
DIAGNOSTIC: Amp×sin sheet / local-density layers — labeled, NOT Gauss.

Cone half-angle: Migkas §4.3 uses "cone of a given radius θ" with θ=75°
(angular separation from center) — confirmed in migkas2020.txt.

Honest caveats:
  - B is mass-proxy AN, not Migkas L_X–T slope.
  - Sine is DIAGNOSTIC, not first-principles.
  - Cone avg is equal-weight (cos-lat solid-angle) idealization, not full
    cluster selection / distance weighting of Migkas.
  - Fig.6 is Δℓ=90° longitude strips (secondary); cone is primary new result.
  - Migkas 2D cone A/B maps (Figs 8–9) are NOT digitised — blocker noted.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

ROOT = Path("/workspace/cosmology")
MAPS = ROOT / "fullsky_AB" / "fullsky_AB_maps.npz"
FIG6 = ROOT / "migkas" / "fig6_digitized.csv"
OUT = Path(__file__).resolve().parent

THETA_DEG = 75.0  # Migkas cone RADIUS (angular separation)
AXES = {
    "sample_280m20": {"l": 280.0, "b": -20.0, "label": "sample (280°, −20°)"},
    "joint_303m27": {"l": 303.0, "b": -27.0, "label": "joint (303°, −27°)"},
}
OMEGA_L = 2.0 * np.pi / 360.0
SINE_VARIANTS = {
    "sheet": {
        "Amp_A": -0.06,
        "Amp_B": +0.06,
        "phi_rad": 0.35,
        "label": "sheet sine DIAGNOSTIC (Amp_A=-0.06, Amp_B=+0.06, φ=0.35)",
    },
    "local_density": {
        "Amp_A": +0.07,
        "Amp_B": -0.07,
        "phi_rad": float(np.deg2rad(212.81767955801106)),
        "label": "local-density sine DIAGNOSTIC (Amp_A=+0.07, Amp_B=-0.07, φ≈212.82°)",
    },
}


def style() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 180,
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "mathtext.fontset": "dejavusans",
        }
    )


def now_local() -> str:
    return datetime.now(ZoneInfo("America/Chicago")).strftime(
        "%Y-%m-%d %I:%M %p %Z (UTC%z)"
    )


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


def rotate_map_to_sky(
    values: np.ndarray,
    Lon: np.ndarray,
    Lat: np.ndarray,
    R: np.ndarray,
    n_lon: int,
    n_lat: int,
) -> np.ndarray:
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


def apply_sine(AF, AN, l_deg, amp_a, amp_b, phi):
    s = np.sin(OMEGA_L * l_deg + phi)
    BA = amp_a * s
    BB = amp_b * s
    return (1.0 - BA) * AF, (1.0 - BB) * AN


def prepare_mollweide(Lon, Lat, Z):
    n_lon = Lon.shape[1]
    i0 = n_lon // 2
    Lon_r = np.roll(Lon, -i0, axis=1)
    Lat_r = np.roll(Lat, -i0, axis=1)
    Z_r = np.roll(Z, -i0, axis=1)
    lon_plot = np.where(Lon_r > np.pi, Lon_r - 2.0 * np.pi, Lon_r)
    return lon_plot, Lat_r, Z_r


def rsq(y_obs: np.ndarray, y_model: np.ndarray) -> float:
    y = np.asarray(y_obs, dtype=float)
    x = np.asarray(y_model, dtype=float)
    if y.size < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1] ** 2)


def boxcar_dl90(y: np.ndarray) -> np.ndarray:
    """Δℓ=90° wrap-around boxcar on a regular 1° longitude grid."""
    n = len(y)
    dl = 360.0 / n
    half = int(round(45.0 / dl))
    kern = np.ones(2 * half + 1)
    kern /= kern.sum()
    return np.convolve(np.pad(y, (half, half), mode="wrap"), kern, mode="valid")


def unit_vectors_from_lb(l_deg: np.ndarray, b_deg: np.ndarray) -> np.ndarray:
    l = np.deg2rad(np.asarray(l_deg, dtype=float))
    b = np.deg2rad(np.asarray(b_deg, dtype=float))
    return np.stack(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], axis=-1
    )


def cone_average_map(
    values: np.ndarray,
    n_pix: np.ndarray,
    w_pix: np.ndarray,
    centers_l: np.ndarray,
    centers_b: np.ndarray,
    theta_deg: float = THETA_DEG,
) -> np.ndarray:
    """
    For each cone center, solid-angle-weighted mean of map pixels with
    angle(n̂, ĉ) < θ. Equal-weight idealization (cos-lat already in w_pix).
    """
    cos_th = np.cos(np.deg2rad(theta_deg))
    n_flat = n_pix.reshape(-1, 3)
    w_flat = w_pix.ravel()
    v_flat = values.ravel()
    # normalize weights once
    # batch centers
    c = unit_vectors_from_lb(centers_l, centers_b)  # (Nc, 3)
    # dots: (Nc, Npix)
    dots = c @ n_flat.T
    out = np.empty(c.shape[0], dtype=float)
    for i in range(c.shape[0]):
        m = dots[i] >= cos_th
        ww = w_flat[m]
        sw = ww.sum()
        if sw <= 0:
            out[i] = np.nan
        else:
            out[i] = np.sum(ww * v_flat[m]) / sw
    return out


def cone_average_strip(
    values, n_pix, w_pix, l_centers, b_fixed, theta_deg=THETA_DEG
):
    b_arr = np.full_like(l_centers, b_fixed, dtype=float)
    return cone_average_map(values, n_pix, w_pix, l_centers, b_arr, theta_deg)


def lat_weighted_lon_profile(values: np.ndarray, lat_1d_rad: np.ndarray) -> np.ndarray:
    """Mean over latitude with cos(lat) weights → A(ℓ) on the sky grid."""
    w = np.cos(lat_1d_rad)  # (n_lat,)
    w = w / w.sum()
    return (values * w[:, None]).sum(axis=0)


def build_center_grid(dl=10.0, db=10.0):
    l = np.arange(0.0, 360.0, dl)
    b = np.arange(-90.0 + db, 90.0, db)  # avoid exact poles
    L, B = np.meshgrid(l, b)
    return L, B


def plot_cone_mollweide(L_c, B_c, maps: dict, path: Path, title: str, disclaimer: str = ""):
    style()
    keys = list(maps.keys())
    nrows = int(np.ceil(len(keys) / 2))
    ncols = 2 if len(keys) > 1 else 1
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(12.5, 3.5 * nrows),
        subplot_kw={"projection": "mollweide"},
    )
    axes = np.atleast_1d(axes).ravel()
    Lon = np.deg2rad(L_c)
    Lat = np.deg2rad(B_c)
    for ax, key in zip(axes, keys):
        Z = maps[key]
        lon_plot, Lat_p, Zp = prepare_mollweide(Lon, Lat, Z)
        pcm = ax.pcolormesh(lon_plot, Lat_p, Zp, shading="auto", cmap="RdBu_r")
        fig.colorbar(pcm, ax=ax, shrink=0.65, pad=0.05)
        ax.set_title(key, fontsize=10)
        ax.grid(True, alpha=0.3)
    for j in range(len(keys), len(axes)):
        axes[j].set_visible(False)
    fig.suptitle(title, fontsize=12, y=1.02)
    if disclaimer:
        fig.text(0.5, -0.02, disclaimer, ha="center", fontsize=8, style="italic",
                 color="#444444", wrap=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_strip_overlay(
    l_model, A_model, B_model,
    migkas,
    path: Path,
    title: str,
    r2a: float, r2b: float,
    disclaimer: str,
    secondary_l=None, secondary_A=None, secondary_B=None, secondary_label=None,
):
    style()
    fig, axes = plt.subplots(2, 1, figsize=(11, 7.2), sharex=True)
    # A
    ax = axes[0]
    ax.errorbar(
        migkas["l_deg"], migkas["A_over_Aall"],
        yerr=[migkas["A_err_low"], migkas["A_err_high"]],
        fmt="o", ms=3.5, color="#222222", ecolor="#888888",
        elinewidth=0.8, capsize=2, label="Migkas Fig.6 A/Aall (digitised)",
        zorder=3,
    )
    ax.plot(l_model, A_model, "-", color="#c0392b", lw=2.0,
            label=f"θ=75° cone avg  (R²={r2a:.3f})")
    if secondary_l is not None:
        ax.plot(secondary_l, secondary_A, "--", color="#e67e22", lw=1.4,
                label=secondary_label or "Δℓ=90° boxcar strip")
    ax.axhline(1.0, color="k", lw=0.6, alpha=0.4)
    ax.set_ylabel(r"$A / A_{\mathrm{all}}$")
    ax.set_ylim(0.7, 1.35)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    ax.errorbar(
        migkas["l_deg"], migkas["B_over_Ball"],
        yerr=[migkas["B_err_low"], migkas["B_err_high"]],
        fmt="o", ms=3.5, color="#222222", ecolor="#888888",
        elinewidth=0.8, capsize=2, label="Migkas Fig.6 B/Ball (digitised)",
        zorder=3,
    )
    ax.plot(l_model, B_model, "-", color="#2471a3", lw=2.0,
            label=f"θ=75° cone avg  (R²={r2b:.3f})")
    if secondary_l is not None:
        ax.plot(secondary_l, secondary_B, "--", color="#5dade2", lw=1.4,
                label=secondary_label or "Δℓ=90° boxcar strip")
    ax.axhline(1.0, color="k", lw=0.6, alpha=0.4)
    ax.set_ylabel(r"$B / B_{\mathrm{all}}$  [model: mass-proxy AN]")
    ax.set_xlabel(r"Galactic longitude $\ell$ [deg]")
    ax.set_ylim(0.7, 1.35)
    ax.set_xlim(0, 360)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.25)

    fig.text(0.5, 0.01, disclaimer, ha="center", fontsize=7.5, style="italic",
             color="#444444")
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_scatter(migkas_A, model_A, migkas_B, model_B, path, title, r2a, r2b):
    style()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, mo, mi, r2, name, c in [
        (axes[0], model_A, migkas_A, r2a, "A", "#c0392b"),
        (axes[1], model_B, migkas_B, r2b, "B", "#2471a3"),
    ]:
        ax.scatter(mo, mi, s=22, c=c, alpha=0.85, edgecolors="k", linewidths=0.3)
        lims = [
            min(mo.min(), mi.min()) - 0.02,
            max(mo.max(), mi.max()) + 0.02,
        ]
        ax.plot(lims, lims, "k--", lw=0.8, alpha=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_aspect("equal", "box")
        ax.set_xlabel(f"cone model {name}")
        ax.set_ylabel(f"Migkas Fig.6 {name}")
        ax.set_title(f"{name}: R² = {r2:.3f}")
        ax.grid(True, alpha=0.25)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def interp_periodic(l_tgt, l_src, y_src):
    order = np.argsort(l_src)
    ls, ys = np.asarray(l_src)[order], np.asarray(y_src)[order]
    ls_p = np.concatenate([ls - 360.0, ls, ls + 360.0])
    ys_p = np.concatenate([ys, ys, ys])
    return np.interp(np.asarray(l_tgt) % 360.0, ls_p, ys_p)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    style()
    print(f"[{now_local()}] loading maps…")
    d = np.load(MAPS)
    Lon = d["lon"]
    Lat = d["lat"]
    lon_1d = d["lon_1d"]
    lat_1d = d["lat_1d"]
    n_lon, n_lat = len(lon_1d), len(lat_1d)
    weights = d["weights"]  # already ∝ cos(lat) dlon dlat-ish

    # Pixel unit vectors in sky/Galactic frame (lon↔ℓ, lat↔b)
    n_pix = np.stack(
        [
            np.cos(Lat) * np.cos(Lon),
            np.cos(Lat) * np.sin(Lon),
            np.sin(Lat),
        ],
        axis=-1,
    )
    l_gal = (np.rad2deg(Lon) % 360.0)

    A_hr_m = d["A_hratio_model"]
    A_inv_m = d["A_invsq_model"]
    B_m = d["B_model"]

    migkas = np.genfromtxt(FIG6, delimiter=",", names=True)
    # structured array fields

    L_c, B_c = build_center_grid(dl=10.0, db=10.0)
    centers_l = L_c.ravel()
    centers_b = B_c.ravel()
    print(f"cone centers: {centers_l.size} (dl=db=10°), θ={THETA_DEG}° radius")

    strip_l = np.arange(0.0, 360.0, 5.0)
    summary = {
        "created_local": now_local(),
        "theta_deg": THETA_DEG,
        "theta_definition": (
            "Migkas cone RADIUS = angular separation from center "
            "(§4.3: 'cone of a given radius θ'; θ=75°)"
        ),
        "weighting": (
            "solid-angle cos(lat) weights from fullsky maps; "
            "equal-weight idealization — NOT Migkas cluster selection / "
            "distance-from-center weighting"
        ),
        "B_caveat": "B is mass-proxy AN = m(r_src)/M_enc mean-renorm, NOT Migkas L_X-T slope",
        "sine_caveat": (
            "sheet Amp×sin is DIAGNOSTIC — not derived from Gauss; "
            "not first-principles gravity"
        ),
        "blocker_2d_cone_maps": (
            "Migkas Figs 8–9 / 21–23 2D cone A/B maps not digitised; "
            "comparison is vs Fig.6 strips + our cone Mollweides only"
        ),
        "axes": {},
        "recommendation": None,
    }

    # Collect R² for winner pick
    scoreboard = []

    for axis_key, axis in AXES.items():
        print(f"\n=== axis {axis_key}: {axis['label']} ===")
        R = rotation_matrix_x_to_sky(axis["l"], axis["b"])
        A_hr = rotate_map_to_sky(A_hr_m, Lon, Lat, R, n_lon, n_lat)
        A_inv = rotate_map_to_sky(A_inv_m, Lon, Lat, R, n_lon, n_lat)
        Bsky = rotate_map_to_sky(B_m, Lon, Lat, R, n_lon, n_lat)

        layers = {
            "physics": {
                "A_hratio": A_hr,
                "A_invsq": A_inv,
                "B": Bsky,
                "kind": "physics",
                "disclaimer": (
                    "PHYSICS-ONLY ST+Gauss; B=mass-proxy AN (not Migkas slope); "
                    "cone avg = cos(lat) equal-weight idealization"
                ),
            }
        }
        for sk, sv in SINE_VARIANTS.items():
            As, Bs = apply_sine(
                A_hr, Bsky, l_gal, sv["Amp_A"], sv["Amp_B"], sv["phi_rad"]
            )
            Ainv_s, _ = apply_sine(
                A_inv, Bsky, l_gal, sv["Amp_A"], sv["Amp_B"], sv["phi_rad"]
            )
            layers[f"sine_{sk}"] = {
                "A_hratio": As,
                "A_invsq": Ainv_s,
                "B": Bs,
                "kind": "sine_diagnostic",
                "disclaimer": (
                    f"DIAGNOSTIC — {sv['label']} — NOT derived from Gauss; "
                    "B=mass-proxy AN; cone avg equal-weight idealization"
                ),
                "sine_label": sv["label"],
            }

        axis_sum = {"label": axis["label"], "l": axis["l"], "b": axis["b"], "layers": {}}

        for layer_key, layer in layers.items():
            print(f"  layer {layer_key}…")
            # --- full cone Mollweide ---
            cone_maps = {}
            minmax = {}
            for mk in ("A_hratio", "A_invsq", "B"):
                avg = cone_average_map(
                    layer[mk], n_pix, weights, centers_l, centers_b, THETA_DEG
                )
                Z = avg.reshape(L_c.shape)
                cone_maps[f"{mk}_cone"] = Z
                minmax[mk] = [float(np.nanmin(Z)), float(np.nanmax(Z))]
                print(f"    {mk} cone min/max = {minmax[mk][0]:.4f} / {minmax[mk][1]:.4f}")

            tag = f"{axis_key}__{layer_key}"
            moll_path = OUT / f"mollweide_cone75_{tag}.png"
            plot_cone_mollweide(
                L_c, B_c, cone_maps, moll_path,
                title=(
                    f"θ=75° cone averages — {axis['label']} — {layer_key}\n"
                    f"{'PHYSICS' if layer['kind']=='physics' else 'DIAGNOSTIC (Amp×sin)'}"
                ),
                disclaimer=layer["disclaimer"],
            )

            # --- strip: cone centers along b = Migkas b ---
            A_strip = cone_average_strip(
                layer["A_hratio"], n_pix, weights, strip_l, axis["b"], THETA_DEG
            )
            B_strip = cone_average_strip(
                layer["B"], n_pix, weights, strip_l, axis["b"], THETA_DEG
            )
            # also b≈0 strip for reference (stored, primary is Migkas b)
            A_strip_b0 = cone_average_strip(
                layer["A_hratio"], n_pix, weights, strip_l, 0.0, THETA_DEG
            )
            B_strip_b0 = cone_average_strip(
                layer["B"], n_pix, weights, strip_l, 0.0, THETA_DEG
            )

            # --- secondary: Δℓ=90° boxcar of lat-weighted lon profile ---
            lon_deg = np.rad2deg(lon_1d) % 360.0
            A_lon = lat_weighted_lon_profile(layer["A_hratio"], lat_1d)
            B_lon = lat_weighted_lon_profile(layer["B"], lat_1d)
            # sort by lon
            order = np.argsort(lon_deg)
            lon_s = lon_deg[order]
            A_box = boxcar_dl90(A_lon[order])
            B_box = boxcar_dl90(B_lon[order])

            # interpolate to Migkas points
            A_at = interp_periodic(migkas["l_deg"], strip_l, A_strip)
            B_at = interp_periodic(migkas["l_deg"], strip_l, B_strip)
            A_box_at = interp_periodic(migkas["l_deg"], lon_s, A_box)
            B_box_at = interp_periodic(migkas["l_deg"], lon_s, B_box)

            r2a = rsq(migkas["A_over_Aall"], A_at)
            r2b = rsq(migkas["B_over_Ball"], B_at)
            r2a_box = rsq(migkas["A_over_Aall"], A_box_at)
            r2b_box = rsq(migkas["B_over_Ball"], B_box_at)
            r2a_b0 = rsq(
                migkas["A_over_Aall"],
                interp_periodic(migkas["l_deg"], strip_l, A_strip_b0),
            )
            r2b_b0 = rsq(
                migkas["B_over_Ball"],
                interp_periodic(migkas["l_deg"], strip_l, B_strip_b0),
            )
            print(
                f"    R² cone@b={axis['b']}: A={r2a:.4f} B={r2b:.4f} | "
                f"boxcar: A={r2a_box:.4f} B={r2b_box:.4f} | "
                f"cone@b=0: A={r2a_b0:.4f} B={r2b_b0:.4f}"
            )

            overlay_path = OUT / f"overlay_fig6_cone75_{tag}.png"
            plot_strip_overlay(
                strip_l, A_strip, B_strip, migkas, overlay_path,
                title=(
                    f"θ=75° cone strip @ b={axis['b']:.0f}° vs Migkas Fig.6 — "
                    f"{axis['label']} — {layer_key}"
                ),
                r2a=r2a, r2b=r2b,
                disclaimer=layer["disclaimer"],
                secondary_l=lon_s, secondary_A=A_box, secondary_B=B_box,
                secondary_label=f"Δℓ=90° boxcar strip (R²A={r2a_box:.3f}, R²B={r2b_box:.3f})",
            )

            scatter_path = OUT / f"scatter_fig6_vs_cone75_{tag}.png"
            plot_scatter(
                migkas["A_over_Aall"], A_at,
                migkas["B_over_Ball"], B_at,
                scatter_path,
                title=f"Migkas Fig.6 vs θ=75° cone — {axis['label']} — {layer_key}",
                r2a=r2a, r2b=r2b,
            )

            # save strip series
            np.savez_compressed(
                OUT / f"strip_series_{tag}.npz",
                strip_l=strip_l,
                A_cone_b_migkas=A_strip,
                B_cone_b_migkas=B_strip,
                A_cone_b0=A_strip_b0,
                B_cone_b0=B_strip_b0,
                lon_box=lon_s,
                A_boxcar90=A_box,
                B_boxcar90=B_box,
                L_centers=L_c,
                B_centers=B_c,
                **{f"cone_{k}": v for k, v in cone_maps.items()},
            )

            layer_sum = {
                "kind": layer["kind"],
                "disclaimer": layer["disclaimer"],
                "minmax_cone": minmax,
                "R2_cone_at_migkas_b": {"A": r2a, "B": r2b, "b_deg": axis["b"]},
                "R2_cone_at_b0": {"A": r2a_b0, "B": r2b_b0},
                "R2_boxcar90_strip": {"A": r2a_box, "B": r2b_box},
                "R2_min_cone": float(min(r2a, r2b)),
                "R2_sum_cone": float(r2a + r2b),
                "files": {
                    "mollweide": str(moll_path),
                    "overlay": str(overlay_path),
                    "scatter": str(scatter_path),
                },
            }
            axis_sum["layers"][layer_key] = layer_sum
            scoreboard.append(
                {
                    "axis": axis_key,
                    "layer": layer_key,
                    "R2_A": r2a,
                    "R2_B": r2b,
                    "R2_min": float(min(r2a, r2b)),
                    "R2_sum": float(r2a + r2b),
                    "kind": layer["kind"],
                }
            )

        summary["axes"][axis_key] = axis_sum

    # Winner picks
    by_sum = sorted(scoreboard, key=lambda r: r["R2_sum"], reverse=True)
    phys = [r for r in scoreboard if r["kind"] == "physics"]
    sine = [r for r in scoreboard if r["kind"] == "sine_diagnostic"]
    best_phys = max(phys, key=lambda r: r["R2_sum"]) if phys else None
    best_sine = max(sine, key=lambda r: r["R2_sum"]) if sine else None
    best_all = by_sum[0]

    summary["scoreboard"] = scoreboard
    summary["best_physics"] = best_phys
    summary["best_sine_diagnostic"] = best_sine
    summary["best_overall"] = best_all
    summary["recommendation"] = (
        f"Overall best Fig.6-strip match: {best_all['layer']} @ {best_all['axis']} "
        f"(R²_A={best_all['R2_A']:.3f}, R²_B={best_all['R2_B']:.3f}). "
        f"Best PHYSICS-ONLY: {best_phys['layer']} @ {best_phys['axis']} "
        f"(R²_A={best_phys['R2_A']:.3f}, R²_B={best_phys['R2_B']:.3f}) — "
        f"expect modest R² (quadrupole-led, two troughs). "
        f"Best SINE DIAGNOSTIC: {best_sine['layer']} @ {best_sine['axis']} "
        f"(R²_A={best_sine['R2_A']:.3f}, R²_B={best_sine['R2_B']:.3f}) — "
        f"Amp×sin is NOT first-principles; for letter use physics cone Mollweides "
        f"as primary geometry result and sine only as labeled diagnostic vs Fig.6. "
        f"Blocker: Migkas 2D θ=75° cone maps not digitised — strip R² is secondary check."
    )

    with open(OUT / "summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    # NOTE.md
    def fmt_row(r):
        return (
            f"| {r['axis']} | {r['layer']} | {r['kind']} | "
            f"{r['R2_A']:.4f} | {r['R2_B']:.4f} | {r['R2_min']:.4f} |"
        )

    rows = "\n".join(fmt_row(r) for r in by_sum)
    # minmax table from physics @ sample
    mm_lines = []
    for ak, axd in summary["axes"].items():
        for lk, ld in axd["layers"].items():
            for mk, (mn, mx) in ld["minmax_cone"].items():
                mm_lines.append(
                    f"| {ak} | {lk} | {mk} | {mn:.4f} | {mx:.4f} |"
                )
    mm_table = "\n".join(mm_lines)

    note = f"""# θ=75° cone scan of full-sky A/B vs Migkas Fig.6

**Date:** {now_local()}

## Cone definition (matched to Migkas)

Migkas+2020 §4.3: *“we consider a cone of a given **radius** θ (we use θ = 45°, 60°,
75° and 90°)”* — θ is the **angular separation** from the cone center (not diameter).
We use **θ = 75°**.

Scanning: cone centers on a lat–lon grid (Δℓ=Δb=10°). For each center ĉ, average
map values over pixels with angle(n̂, ĉ) < 75°, with **cos(lat) solid-angle weights**
from the parent full-sky maps. Equal-weight idealization — **not** Migkas cluster
selection or distance-from-center statistical weights.

## Honest labels

| Item | Meaning here |
|---|---|
| **A_hratio** | Physics AF = H_obs/H mean-renorm (primary) |
| **A_invsq** | Physics AF = 1/(H_obs/H)² mean-renorm (secondary) |
| **B** | Mass-proxy AN = m(r_src)/M_enc mean-renorm — **NOT** Migkas L_X–T slope B |
| **sine_*** | Amp×sin sheet / local-density **DIAGNOSTIC** — **not** derived from Gauss |
| **cone avg** | Idealized solid-angle mean — not full cluster pipeline |

## Axes scanned

- Sample preferred: (ℓ, b) = (280°, −20°)
- Joint 3-sample: (ℓ, b) = (303°, −27°)

Physics maps re-rotated from model frame for each axis. Sine diagnostic applied
as sky function of Galactic ℓ after rotation.

## R² vs digitised Migkas Fig.6 (θ=75° cone strip @ Migkas b)

Primary comparison: A_cone(ℓ), B_cone(ℓ) for centers along b = Migkas b, vs Fig.6.
Secondary: Δℓ=90° boxcar of latitude-weighted full-sky longitude profile (known strip match).

| axis | layer | kind | R²_A | R²_B | R²_min |
|---|---|---|---|---|---|
{rows}

### Best PHYSICS-ONLY
`{best_phys['layer']}` @ `{best_phys['axis']}` — R²_A={best_phys['R2_A']:.4f}, R²_B={best_phys['R2_B']:.4f}

### Best SINE DIAGNOSTIC
`{best_sine['layer']}` @ `{best_sine['axis']}` — R²_A={best_sine['R2_A']:.4f}, R²_B={best_sine['R2_B']:.4f}

### Overall (includes diagnostic)
`{best_all['layer']}` @ `{best_all['axis']}` — R²_A={best_all['R2_A']:.4f}, R²_B={best_all['R2_B']:.4f}

## Cone A/B min/max (over cone-center sky)

| axis | layer | map | min | max |
|---|---|---|---|---|
{mm_table}

## Blocker

Migkas **2D** θ=75° cone A/B significance maps (Figs. 8–9, 21–23) are **not digitised**.
This deliverable compares to **Fig.6 longitude strips** + produces our own cone Mollweides.
Digitising those 2D maps (or re-running on the public cluster table) is the next step for
a true apples-to-apples cone↔cone comparison.

## Recommendation for letter

{summary['recommendation']}

## Regenerable

```bash
/workspace/cosmology/.venv/bin/python /workspace/cosmology/cone_scan_75/run_cone_scan.py
```

Inputs: `/workspace/cosmology/fullsky_AB/fullsky_AB_maps.npz`,
`/workspace/cosmology/migkas/fig6_digitized.csv`.
"""
    (OUT / "NOTE.md").write_text(note)
    print("\n" + summary["recommendation"])
    print(f"\nWrote outputs under {OUT}")
    return summary


if __name__ == "__main__":
    main()

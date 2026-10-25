#!/usr/bin/env python3
"""
Full-sky first-principles A(n̂) and B(n̂) from Chris ST+Gauss two-scale geometry.

PRIMARY: physics H_obs only (no Amp×sin) — first-principles product.
DIAGNOSTIC: sheet sine rotation (Amp×sin) applied as a labeled sky-function layer;
            NOT derived from Gauss / not first-principles gravity.

Binding baseline: twoscale f_u=0.94, D/L=5.9, κ≈0.231, ρ(D)=0.3 ρ_crit.
AF: hratio (primary) + invsq (secondary); B: AN = m(r_src)/M_enc mean-renorm.
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
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "density_local_match"))

from forward_model import make_latlon_grid, fit_dipole_quadrupole  # noqa: E402
from run_density_local_match import (  # noqa: E402
    params_twoscale_local_match,
    peculiar_u,
    m_enc_at_r,
    D,
    H,
    ALPHA,
    EPS,
    OMEGA_M,
)

OUT = Path(__file__).resolve().parent

# Binding baseline
X_DL = 5.9
F_U = 0.94
N_LON = 360
N_LAT = 181

# Migkas-ish preferred axis for sky-frame display (Galactic l,b)
MIGKAS_L_DEG = 280.0
MIGKAS_B_DEG = -20.0

# Sheet sine-rotation diagnostic conventions (NOT first-principles)
# Chris sheet defaults + local-density balanced fit
SINE_VARIANTS = {
    "sheet": {
        "Amp_A": -0.06,
        "Amp_B": +0.06,  # opposite Amp for B
        "phi_rad": 0.35,
        "label": "sheet (Amp_A=-0.06, Amp_B=+0.06, φ=0.35 rad)",
    },
    "local_density": {
        "Amp_A": +0.07,
        "Amp_B": -0.07,  # opposite Amp for B (sheet convention)
        "phi_rad": float(np.deg2rad(212.81767955801106)),
        "label": "local-density fit (Amp_A=+0.07, Amp_B=-0.07, φ≈212.82°)",
    },
}

OMEGA_L = 2.0 * np.pi / 360.0  # sheet: sin(2π ℓ/360 + φ)


def style() -> None:
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.dpi": 200,
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "mathtext.fontset": "dejavusans",
        }
    )


def H_obs_sphere(pmeta: dict, eps: float, n_hat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """3D extension of H_obs_circle. CoE at origin; observer at +D x̂; D̂ = +x."""
    d = eps * D
    D_vec = np.array([D, 0.0, 0.0])
    r_vec = D_vec + d * n_hat
    u_r = peculiar_u(r_vec, pmeta, alpha=ALPHA)
    u_D = peculiar_u(D_vec.reshape(1, 3), pmeta, alpha=ALPHA)[0]
    delta_v = np.sum((u_r - u_D) * n_hat, axis=-1)
    H_obs = H + delta_v / d
    m_src = m_enc_at_r(np.linalg.norm(r_vec, axis=-1), pmeta)
    return H_obs, m_src


def rotation_matrix_x_to_sky(l_deg: float, b_deg: float) -> np.ndarray:
    """R such that R @ e_x = n_gal(l,b); complete triad with Galactic-north preference."""
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
    """Sky value at n_sky = model value at R^T n_sky (bilinear resample)."""
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


def galactic_lon_of_sky_grid(Lon: np.ndarray, Lat: np.ndarray) -> np.ndarray:
    """
    Galactic longitude ℓ (deg) for each sky-frame pixel.
    Sky-frame Lon/Lat use the same Cartesian convention as n_gal:
    n = (cos lat cos lon, cos lat sin lon, sin lat) with lon↔ℓ, lat↔b.
    So ℓ_deg = lon_deg (mod 360) on the sky-frame grid.
    """
    return (np.rad2deg(Lon) % 360.0)


def apply_sine_rotation(
    AF: np.ndarray,
    AN: np.ndarray,
    l_deg: np.ndarray,
    amp_a: float,
    amp_b: float,
    phi: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Sheet convention as sky function of Galactic longitude:
      BA = Amp_A * sin(2π ℓ/360 + φ)
      BB = Amp_B * sin(2π ℓ/360 + φ)   (Amp_B = −Amp_A typically)
      Model A = (1 − BA) * AF
      Model B = (1 − BB) * AN
    Returns (A_mod, B_mod, BA).
    """
    s = np.sin(OMEGA_L * l_deg + phi)
    BA = amp_a * s
    BB = amp_b * s
    return (1.0 - BA) * AF, (1.0 - BB) * AN, BA


def prepare_mollweide_coords(Lon: np.ndarray, Lat: np.ndarray, Z: np.ndarray):
    """Roll so longitude is monotonic in [-π, π] for matplotlib pcolormesh."""
    n_lon = Lon.shape[1]
    i0 = n_lon // 2  # Lon[:, i0] ≈ π
    Lon_r = np.roll(Lon, -i0, axis=1)
    Lat_r = np.roll(Lat, -i0, axis=1)
    Z_r = np.roll(Z, -i0, axis=1)
    lon_plot = np.where(Lon_r > np.pi, Lon_r - 2.0 * np.pi, Lon_r)
    return lon_plot, Lat_r, Z_r


def plot_mollweide_2x2(
    Lon: np.ndarray,
    Lat: np.ndarray,
    maps: dict[str, np.ndarray],
    titles: dict[str, str],
    path: Path,
    frame_label: str,
    galactic_labels: bool = False,
    subtitle_extra: str = "",
) -> None:
    style()
    keys = list(maps.keys())
    nrows = int(np.ceil(len(keys) / 2))
    ncols = 2 if len(keys) > 1 else 1
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(12.5, 3.6 * nrows),
        subplot_kw={"projection": "mollweide"},
    )
    axes = np.atleast_1d(axes).ravel()
    for ax, key in zip(axes, keys):
        lon_plot, Lat_p, Z = prepare_mollweide_coords(Lon, Lat, maps[key])
        pcm = ax.pcolormesh(lon_plot, Lat_p, Z, shading="auto", cmap="RdBu_r")
        if "frac" in key or key.startswith("BA") or "modulation" in key:
            vmax = max(abs(float(Z.min())), abs(float(Z.max())), 1e-6)
            pcm.set_clim(-vmax, vmax)
        else:
            d = max(abs(float(Z.min()) - 1.0), abs(float(Z.max()) - 1.0), 1e-6)
            pcm.set_clim(1.0 - d, 1.0 + d)
        cb = fig.colorbar(
            pcm, ax=ax, orientation="horizontal", pad=0.08, fraction=0.05, aspect=30
        )
        cb.ax.tick_params(labelsize=8)
        ax.set_title(titles.get(key, key), pad=10)
        ax.grid(True, linestyle=":", alpha=0.5)
        if galactic_labels:
            ax.set_xlabel(r"Galactic $\ell$")
        else:
            ax.set_xlabel(r"model lon (from $+x$)")
    for ax in axes[len(keys):]:
        ax.set_visible(False)
    fig.suptitle(
        f"{frame_label}\n"
        rf"twoscale $f_u={F_U}$, $D/L={X_DL}$, $\alpha={ALPHA}$, $d/D={EPS}$"
        + (f"\n{subtitle_extra}" if subtitle_extra else ""),
        fontsize=11,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_single_mollweide(
    Lon, Lat, Z, title, path, center=1.0, galactic=False
) -> None:
    style()
    lon_plot, Lat_p, Zp = prepare_mollweide_coords(Lon, Lat, Z)
    fig = plt.figure(figsize=(10.5, 5.2))
    ax = fig.add_subplot(111, projection="mollweide")
    pcm = ax.pcolormesh(lon_plot, Lat_p, Zp, shading="auto", cmap="RdBu_r")
    Z = Zp  # for clim below
    if center == 0.0:
        vmax = max(abs(float(Z.min())), abs(float(Z.max())), 1e-6)
        pcm.set_clim(-vmax, vmax)
    else:
        d = max(abs(float(Z.min()) - center), abs(float(Z.max()) - center), 1e-6)
        pcm.set_clim(center - d, center + d)
    fig.colorbar(pcm, ax=ax, orientation="horizontal", pad=0.08, fraction=0.06)
    ax.set_title(title, pad=12)
    ax.grid(True, linestyle=":", alpha=0.5)
    if galactic:
        ax.set_xlabel(r"Galactic longitude $\ell$")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_note(meta: dict, path: Path) -> None:
    tz = ZoneInfo("America/Chicago")
    now = datetime.now(tz).strftime("%Y-%m-%d %I:%M %p %Z")
    mp = meta["multipoles_frac"]
    text = f"""# Full-sky A(n̂) and B(n̂) — first principles + sheet sine diagnostic

**Date:** {now}

## What A and B mean here (PRIMARY — first principles)

These maps are **physics-only** angular series from Chris’s ST+Gauss two-scale geometry
(binding baseline: `twoscale_fu0.94`, D/L=5.9, f_u=0.94, κ≈0.231, ρ(D)=0.3 ρ_crit).

| Symbol | Definition | Mean-renorm |
|---|---|---|
| **A_hratio** (primary AF) | A_raw = H_obs/H | A / ⟨A⟩ |
| **A_invsq** (secondary) | A_raw = 1/(H_obs/H)² (H_obs>0) | A / ⟨A⟩ |
| **B** | AN = m(r_src)/M_enc | B / ⟨B⟩ |
| **frac** | H_obs/H − 1 | (not renormed) |

**Honest:** B here is the **mass-proxy AN** used in `density_local_match`, **not** the Migkas
L_X–T slope parameter B from the observational fit. Do not conflate them.

Primary Mollweides (`mollweide_AB_*_frame.png`, `mollweide_sky_*.png`) use **raw physics
H_obs and m(r) only** — Amp×sin is **not** mixed into those maps.

## Sheet sine rotation — DIAGNOSTIC LAYER (retained per Chris)

Chris wants the unexplained sheet Amp×sin phase rotation **kept** as a labeled
diagnostic. It is **not** derived from Gauss / ST and is **not** first-principles gravity.

Sheet convention (applied as a sky function of Galactic longitude ℓ of each pixel):

```
BA = Amp_A · sin(2π ℓ/360 + φ)
BB = Amp_B · sin(2π ℓ/360 + φ)   # Amp_B = −Amp_A (opposite)
Model A = (1 − BA) · AF
Model B = (1 − BB) · AN
```

Two diagnostic variants (sky-frame AF=A_hratio, AN=B):

1. **sheet defaults:** Amp_A=−0.06, Amp_B=+0.06, φ=0.35 rad
2. **local-density balanced fit:** Amp_A=+0.07, Amp_B=−0.07, φ≈212.82°

Outputs clearly labeled:
`mollweide_sine_*_DIAGNOSTIC.png` — **"sheet sine rotation diagnostic — not derived from Gauss"**.

Do **not** present sine-rotated maps as first-principles gravity. They are an
axis/modulation diagnostic Chris wants retained for Migkas Fig.6-style comparisons.

## Quadrupole-led expectation (physics frac)

At d/D = ε = 0.1, the tidal limit of ST+Gauss gives Ψ/H ≈ ακ(3μ²−1) = 2ακ P₂(μ):
**pure quadrupole, vanishing dipole at leading order.** Finite ε generates a small
dipole. Multipole moments on `frac` (model frame):

- |A_dipole| = {mp['A']:.4e}
- Q_P2 = {mp['Q_P2_moment']:.4e}
- |C1| ≪ |C2|: |A| / |Q_P2| = {abs(mp['A'])/max(abs(mp['Q_P2_moment']),1e-30):.4e}

## Frames and rotation

- **Model frame:** CoE at origin, observer at +D x̂, preferred axis D̂ = +x.
  Lon measured from +x in the xy-plane; lat = elevation.
- **Sky frame (display):** rotate so +x → Migkas-ish preferred direction
  (ℓ, b) ≈ ({MIGKAS_L_DEG}°, {MIGKAS_B_DEG}°).
  Rotation R maps e_x → n_gal(ℓ,b); sky value at n_sky = model value at Rᵀ n_sky.
  Documented in `metadata.json` (`rotation_matrix`, `migkas_lb_deg`).
- Sine diagnostic uses **Galactic ℓ of each sky-frame pixel** (lon on the sky grid).

## Baseline parameters

- profile: twoscale, D/L={X_DL}, f_u={F_U}
- κ = {meta['pmeta']['kappa']:.6f}, κ_core = {meta['pmeta']['kappa_core']:.6f}
- α = {ALPHA}, d/D = {EPS}, Ω_m = {OMEGA_M}
- grid: n_lon={N_LON}, n_lat={N_LAT}

## Next step

θ = 75° cone average (directional averages of A/B over Migkas-style cones)
for direct comparison to Fig.6-style longitude series. Cone averages can be run
on physics-only maps and, separately, on the labeled sine-diagnostic maps.

## Regenerable

```bash
/workspace/cosmology/.venv/bin/python /workspace/cosmology/fullsky_AB/run_fullsky_AB.py
```
"""
    path.write_text(text)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pmeta = params_twoscale_local_match(X_DL, F_U)
    assert abs(pmeta["rho_local_over_rcrit"] - OMEGA_M) < 1e-8
    print("pmeta:", json.dumps({k: pmeta[k] for k in (
        "profile", "x", "f_u", "kappa", "kappa_core", "M_enc",
        "rho_local_over_rcrit", "rho_mean_over_rcrit",
    )}, indent=2), flush=True)

    n_hat, lon_flat, lat_flat, w = make_latlon_grid(N_LON, N_LAT)
    Lon = lon_flat.reshape(N_LAT, N_LON)
    Lat = lat_flat.reshape(N_LAT, N_LON)

    print(f"Computing H_obs on {N_LON}×{N_LAT} grid …", flush=True)
    H_obs, m_src = H_obs_sphere(pmeta, EPS, n_hat)
    H_safe = np.maximum(H_obs, 1e-6 * H)
    frac = H_obs / H - 1.0

    A_hratio_raw = H_safe / H
    A_invsq_raw = 1.0 / (H_safe / H) ** 2
    B_raw = m_src / pmeta["M_enc"]

    def mean_renorm(x):
        m = float(np.average(x, weights=w))
        return x / m, m

    A_hratio, mean_Ah = mean_renorm(A_hratio_raw)
    A_invsq, mean_Ai = mean_renorm(A_invsq_raw)
    B, mean_B = mean_renorm(B_raw)

    fit_frac = fit_dipole_quadrupole(frac, n_hat, w)
    fit_Ah = fit_dipole_quadrupole(A_hratio - 1.0, n_hat, w)
    fit_B = fit_dipole_quadrupole(B - 1.0, n_hat, w)

    print("=== Multipoles on frac (model frame) ===", flush=True)
    print(f"  monopole={fit_frac['monopole']:.6e}", flush=True)
    print(f"  A_dipole={fit_frac['A']:.6e}  A_vec={fit_frac['A_vec']}", flush=True)
    print(f"  A_x_moment={fit_frac['A_x_moment']:.6e}", flush=True)
    print(f"  Q_P2_moment={fit_frac['Q_P2_moment']:.6e}", flush=True)
    print(
        f"  |A|/|Q_P2|={abs(fit_frac['A'])/max(abs(fit_frac['Q_P2_moment']),1e-30):.4e}  "
        f"(expect ≪1 at d/D={EPS})",
        flush=True,
    )
    print(f"  frac min/max/ptp = {frac.min():.6f}/{frac.max():.6f}/{np.ptp(frac):.6f}", flush=True)
    print(f"  A_hratio min/max = {A_hratio.min():.6f}/{A_hratio.max():.6f}", flush=True)
    print(f"  A_invsq  min/max = {A_invsq.min():.6f}/{A_invsq.max():.6f}", flush=True)
    print(f"  B        min/max = {B.min():.6f}/{B.max():.6f}", flush=True)

    maps_model = {
        "A_hratio": A_hratio.reshape(N_LAT, N_LON),
        "A_invsq": A_invsq.reshape(N_LAT, N_LON),
        "B": B.reshape(N_LAT, N_LON),
        "frac": frac.reshape(N_LAT, N_LON),
        "H_obs": H_obs.reshape(N_LAT, N_LON),
    }

    R = rotation_matrix_x_to_sky(MIGKAS_L_DEG, MIGKAS_B_DEG)
    n_check = R @ np.array([1.0, 0.0, 0.0])
    l_chk = np.rad2deg(np.arctan2(n_check[1], n_check[0])) % 360.0
    b_chk = np.rad2deg(np.arcsin(np.clip(n_check[2], -1, 1)))
    print(
        f"Rotation check: +x → (l,b)=({l_chk:.2f}°, {b_chk:.2f}°) "
        f"[target {MIGKAS_L_DEG}, {MIGKAS_B_DEG}]",
        flush=True,
    )

    maps_sky = {
        k: rotate_map_to_sky(maps_model[k].ravel(), Lon, Lat, R, N_LON, N_LAT)
        for k in ("A_hratio", "A_invsq", "B", "frac")
    }

    titles_phys = {
        "A_hratio": r"$A_{\mathrm{hratio}}=(H_{\mathrm{obs}}/H)/\mathrm{mean}$ [physics]",
        "A_invsq": r"$A_{\mathrm{invsq}}\propto 1/(H_{\mathrm{obs}}/H)^2$ [physics]",
        "B": r"$B_{\mathrm{AN}}=m(r_{\mathrm{src}})/M_{\mathrm{enc}}$ [physics]",
        "frac": r"$\mathrm{frac}=H_{\mathrm{obs}}/H-1$ [physics]",
    }

    # --- PRIMARY: physics-only Mollweides ---
    plot_mollweide_2x2(
        Lon, Lat, {k: maps_model[k] for k in titles_phys}, titles_phys,
        OUT / "mollweide_AB_model_frame.png",
        frame_label="PHYSICS-ONLY — model frame (CoE / D̂ on +x)",
        galactic_labels=False,
        subtitle_extra=r"first principles; no Amp$\times$sin",
    )
    plot_mollweide_2x2(
        Lon, Lat, {k: maps_sky[k] for k in titles_phys}, titles_phys,
        OUT / "mollweide_AB_sky_frame.png",
        frame_label=rf"PHYSICS-ONLY — sky (pref. axis $\to$ $\ell={MIGKAS_L_DEG:.0f}^\circ$, $b={MIGKAS_B_DEG:.0f}^\circ$)",
        galactic_labels=True,
        subtitle_extra=r"first principles; no Amp$\times$sin",
    )
    for key, center in [("A_hratio", 1.0), ("A_invsq", 1.0), ("B", 1.0), ("frac", 0.0)]:
        plot_single_mollweide(
            Lon, Lat, maps_sky[key],
            titles_phys[key]
            + rf" — sky ($\ell={MIGKAS_L_DEG:.0f}^\circ$, $b={MIGKAS_B_DEG:.0f}^\circ$)",
            OUT / f"mollweide_sky_{key}.png",
            center=center,
            galactic=True,
        )

    # --- DIAGNOSTIC: sheet sine rotation on sky-frame AF/AN ---
    l_gal = galactic_lon_of_sky_grid(Lon, Lat)
    sine_maps_saved = {}
    sine_meta = {}
    AF_sky = maps_sky["A_hratio"]
    AN_sky = maps_sky["B"]
    AF_inv_sky = maps_sky["A_invsq"]

    for vname, vcfg in SINE_VARIANTS.items():
        print(f"Sine diagnostic: {vname} — {vcfg['label']}", flush=True)
        A_mod, B_mod, BA = apply_sine_rotation(
            AF_sky, AN_sky, l_gal,
            vcfg["Amp_A"], vcfg["Amp_B"], vcfg["phi_rad"],
        )
        A_inv_mod, _, _ = apply_sine_rotation(
            AF_inv_sky, AN_sky, l_gal,
            vcfg["Amp_A"], vcfg["Amp_B"], vcfg["phi_rad"],
        )
        maps_sine = {
            "A_hratio_sine": A_mod,
            "A_invsq_sine": A_inv_mod,
            "B_sine": B_mod,
            "BA_modulation": BA,
        }
        titles_sine = {
            "A_hratio_sine": r"$(1-BA)\,A_{\mathrm{hratio}}$ [DIAGNOSTIC]",
            "A_invsq_sine": r"$(1-BA)\,A_{\mathrm{invsq}}$ [DIAGNOSTIC]",
            "B_sine": r"$(1-BB)\,B_{\mathrm{AN}}$ [DIAGNOSTIC]",
            "BA_modulation": r"$BA=\mathrm{Amp}_A\sin(2\pi\ell/360+\phi)$",
        }
        plot_mollweide_2x2(
            Lon, Lat, maps_sine, titles_sine,
            OUT / f"mollweide_sine_{vname}_DIAGNOSTIC.png",
            frame_label=(
                "DIAGNOSTIC — sheet sine rotation "
                "(NOT derived from Gauss / NOT first-principles)"
            ),
            galactic_labels=True,
            subtitle_extra=vcfg["label"],
        )
        plot_single_mollweide(
            Lon, Lat, A_mod,
            rf"DIAGNOSTIC $(1-BA)A_{{\mathrm{{hratio}}}}$ — {vcfg['label']}",
            OUT / f"mollweide_sine_{vname}_A_hratio_DIAGNOSTIC.png",
            center=1.0, galactic=True,
        )
        plot_single_mollweide(
            Lon, Lat, B_mod,
            rf"DIAGNOSTIC $(1-BB)B_{{\mathrm{{AN}}}}$ — {vcfg['label']}",
            OUT / f"mollweide_sine_{vname}_B_DIAGNOSTIC.png",
            center=1.0, galactic=True,
        )
        sine_maps_saved[vname] = {
            "A_hratio_sine": A_mod,
            "A_invsq_sine": A_inv_mod,
            "B_sine": B_mod,
            "BA": BA,
        }
        sine_meta[vname] = {
            **{k: vcfg[k] for k in ("Amp_A", "Amp_B", "phi_rad", "label")},
            "phi_deg": float(np.rad2deg(vcfg["phi_rad"])),
            "minmax_A_hratio_sine": [float(A_mod.min()), float(A_mod.max())],
            "minmax_B_sine": [float(B_mod.min()), float(B_mod.max())],
            "disclaimer": (
                "sheet sine rotation diagnostic — not derived from Gauss; "
                "not first-principles gravity"
            ),
        }

    tz = ZoneInfo("America/Chicago")
    now = datetime.now(tz)
    meta = {
        "created_local": now.strftime("%Y-%m-%d %I:%M %p %Z (UTC%z)"),
        "baseline": {
            "profile": "twoscale_fu0.94",
            "D_over_L": X_DL,
            "f_u": F_U,
            "alpha": ALPHA,
            "d_over_D": EPS,
            "Omega_m": OMEGA_M,
            "af_modes": ["hratio", "invsq"],
            "B_proxy": "AN = m(r_src)/M_enc mean-renorm",
            "primary_product": "physics H_obs / m(r) only — Amp×sin NOT in primary maps",
            "sine_diagnostic_retained": True,
            "sine_disclaimer": (
                "sheet sine rotation diagnostic — not derived from Gauss; "
                "not first-principles gravity"
            ),
        },
        "pmeta": pmeta,
        "grid": {"n_lon": N_LON, "n_lat": N_LAT},
        "migkas_lb_deg": [MIGKAS_L_DEG, MIGKAS_B_DEG],
        "rotation_matrix": R.tolist(),
        "rotation_doc": (
            "R maps model +x (D̂) to Galactic n(l,b)=(cos b cos l, cos b sin l, sin b). "
            "Sky-frame map value at n_sky equals model-frame value at R^T n_sky. "
            "Sine diagnostic uses Galactic ℓ = sky-frame lon (deg) of each pixel."
        ),
        "means_raw": {
            "A_hratio": mean_Ah,
            "A_invsq": mean_Ai,
            "B": mean_B,
        },
        "minmax": {
            "A_hratio": [float(A_hratio.min()), float(A_hratio.max())],
            "A_invsq": [float(A_invsq.min()), float(A_invsq.max())],
            "B": [float(B.min()), float(B.max())],
            "frac": [float(frac.min()), float(frac.max())],
        },
        "multipoles_frac": fit_frac,
        "multipoles_A_hratio_minus1": fit_Ah,
        "multipoles_B_minus1": fit_B,
        "sine_diagnostic": sine_meta,
    }
    (OUT / "metadata.json").write_text(json.dumps(meta, indent=2))

    save_kw = dict(
        lon=Lon,
        lat=Lat,
        lon_1d=np.linspace(0.0, 2.0 * np.pi, N_LON, endpoint=False),
        lat_1d=np.linspace(-0.5 * np.pi + 1e-3, 0.5 * np.pi - 1e-3, N_LAT),
        A_hratio_model=maps_model["A_hratio"],
        A_invsq_model=maps_model["A_invsq"],
        B_model=maps_model["B"],
        frac_model=maps_model["frac"],
        H_obs_model=maps_model["H_obs"],
        A_hratio_sky=maps_sky["A_hratio"],
        A_invsq_sky=maps_sky["A_invsq"],
        B_sky=maps_sky["B"],
        frac_sky=maps_sky["frac"],
        l_gal_sky_deg=l_gal,
        weights=w.reshape(N_LAT, N_LON),
        rotation_matrix=R,
        n_hat=n_hat.reshape(N_LAT, N_LON, 3),
    )
    for vname, sm in sine_maps_saved.items():
        for k, arr in sm.items():
            save_kw[f"{k}_{vname}"] = arr
    np.savez_compressed(OUT / "fullsky_AB_maps.npz", **save_kw)

    write_note(meta, OUT / "NOTE.md")

    summary = {
        "minmax": meta["minmax"],
        "Q_P2_frac": fit_frac["Q_P2_moment"],
        "A_dipole_frac": fit_frac["A"],
        "A_x_moment_frac": fit_frac["A_x_moment"],
        "ratio_A_over_Q": abs(fit_frac["A"]) / max(abs(fit_frac["Q_P2_moment"]), 1e-30),
        "sine_diagnostic_variants": list(SINE_VARIANTS.keys()),
        "sine_disclaimer": meta["baseline"]["sine_disclaimer"],
        "outputs": {
            "npz": str(OUT / "fullsky_AB_maps.npz"),
            "metadata": str(OUT / "metadata.json"),
            "note": str(OUT / "NOTE.md"),
            "mollweide_model_physics": str(OUT / "mollweide_AB_model_frame.png"),
            "mollweide_sky_physics": str(OUT / "mollweide_AB_sky_frame.png"),
            "mollweide_sine_sheet": str(OUT / "mollweide_sine_sheet_DIAGNOSTIC.png"),
            "mollweide_sine_local_density": str(
                OUT / "mollweide_sine_local_density_DIAGNOSTIC.png"
            ),
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print("Wrote outputs under", OUT, flush=True)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()

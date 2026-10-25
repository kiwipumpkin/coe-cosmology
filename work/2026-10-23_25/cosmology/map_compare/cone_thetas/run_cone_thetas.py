#!/usr/bin/env python3
"""4×4 grid of cone-averaged A/A_i Mollweides at Migkas cone sizes θ=90,75,60,45."""
from __future__ import annotations

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
sys.path.insert(0, str(ROOT / "map_compare"))

from run_levelset_maps import (  # noqa: E402
    wrap_lon_deg,
    lon_disp_from_ell,
    celestial_equator_lb,
    style,
    VMIN,
    VMAX,
)

OUT = Path(__file__).resolve().parent
FULLSKY = ROOT / "fullsky_AB" / "fullsky_AB_maps.npz"
CONE_REF = ROOT / "map_compare" / "cone75_maps.npz"

THETAS = [90.0, 75.0, 60.0, 45.0]
COL_KEYS = [
    ("phys_invsq", "A_invsq_sky", "Physics invsq"),
    ("sheet_invsq", "A_invsq_sine_sheet", "Sheet invsq"),
    ("phys_hratio", "A_hratio_sky", "Physics hratio"),
    ("sheet_hratio", "A_hratio_sine_sheet", "Sheet hratio"),
]


def now_local() -> str:
    return datetime.now(ZoneInfo("America/Chicago")).strftime(
        "%Y-%m-%d %I:%M %p %Z (UTC%z)"
    )


def build_cone_geometry(L_centers, B_centers, n_hat_pix, weights, theta_deg: float):
    Lc, Bc = np.meshgrid(L_centers, B_centers)
    l = np.deg2rad(Lc.ravel())
    b = np.deg2rad(Bc.ravel())
    c = np.stack(
        [np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], axis=-1
    )
    n_flat = n_hat_pix.reshape(-1, 3)
    w_flat = weights.ravel()
    cos_th = np.cos(np.deg2rad(theta_deg))
    return c, n_flat, w_flat, cos_th, Lc.shape


def cone_average(values, geom):
    c, n_flat, w_flat, cos_th, shape = geom
    v_flat = values.ravel()
    out = np.empty(c.shape[0], dtype=float)
    chunk = 200
    for i0 in range(0, c.shape[0], chunk):
        i1 = min(i0 + chunk, c.shape[0])
        dots = c[i0:i1] @ n_flat.T
        for j in range(i1 - i0):
            m = dots[j] >= cos_th
            ww = w_flat[m]
            sw = ww.sum()
            out[i0 + j] = np.sum(ww * v_flat[m]) / sw if sw > 0 else np.nan
    Z = out.reshape(shape)
    return Z / np.nanmean(Z)


def prep_mesh(L_c, B_c, Z):
    Lon_c, Lat_c = np.meshgrid(L_c, B_c)
    lon_disp = np.deg2rad(lon_disp_from_ell(Lon_c))
    lat = np.deg2rad(Lat_c)
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
        ax.plot(ld[s:e], bd[s:e], "k-", lw=0.9, alpha=0.8, zorder=5)


def plot_panel(ax, L_c, B_c, Z, title=None, xticklabels=True):
    Lon_e, Lat_e, Z_s = prep_mesh(L_c, B_c, Z)
    pcm = ax.pcolormesh(
        Lon_e,
        Lat_e,
        Z_s,
        shading="flat",
        cmap="magma",
        norm=Normalize(VMIN, VMAX),
        rasterized=True,
    )
    ax.grid(True, linestyle=":", alpha=0.4, color="white", lw=0.45)
    draw_equator(ax)
    ax.set_xticks(np.deg2rad([-90, 0, 90]))
    if xticklabels:
        ax.set_xticklabels([r"$90^\circ$", r"$0^\circ$", r"$270^\circ$"], fontsize=8)
    else:
        ax.set_xticklabels([])
    ax.set_yticks([])
    if title:
        ax.set_title(title, fontsize=11, pad=6)
    return pcm


def plot_card(L_c, B_c, Z, title, path):
    style()
    fig = plt.figure(figsize=(10, 6.5), dpi=120)
    ax = fig.add_subplot(111, projection="mollweide")
    pcm = plot_panel(ax, L_c, B_c, Z, title=title)
    cbar = fig.colorbar(pcm, ax=ax, orientation="vertical", pad=0.02, fraction=0.035, aspect=28)
    ticks = np.linspace(VMIN, VMAX, 8)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t:.2f}" for t in ticks])
    cbar.set_label(r"$A/A_i$", fontsize=11)
    fig.subplots_adjust(left=0.04, right=0.92, bottom=0.06, top=0.88)
    fig.savefig(path, dpi=120, facecolor="white")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cards_dir = OUT / "cards"
    cards_dir.mkdir(exist_ok=True)
    t0 = time.time()
    print(f"[{now_local()}] cone θ grid → {OUT}", flush=True)

    ref = np.load(CONE_REF)
    L_c = np.asarray(ref["L_centers"], dtype=float)
    B_c = np.asarray(ref["B_centers"], dtype=float)
    print(f"grid L={L_c.shape} B={B_c.shape}", flush=True)

    fs = np.load(FULLSKY)
    n_hat = fs["n_hat"]
    weights = fs["weights"]
    sky_maps = {k: fs[src] for _, src, _ in COL_KEYS for k in [src]}
    # ensure we have the four
    for _, src, _ in COL_KEYS:
        assert src in fs.files, src

    maps = {}  # (theta, col_key) -> Z
    ranges = {}

    for theta in THETAS:
        print(f"\n=== θ={theta:.0f}° ===", flush=True)
        geom = build_cone_geometry(L_c, B_c, n_hat, weights, theta)
        for col_key, src, _ in COL_KEYS:
            t1 = time.time()
            Z = cone_average(fs[src], geom)
            maps[(theta, col_key)] = Z
            mn, mx = float(np.nanmin(Z)), float(np.nanmax(Z))
            ranges[(theta, col_key)] = (mn, mx)
            print(
                f"  {col_key:14s} [{mn:.4f}, {mx:.4f}]  ({time.time()-t1:.1f}s)",
                flush=True,
            )
            # optional individual cards
            title = rf"θ={theta:.0f}° cone $A/A_i$ — {col_key.replace('_', ' ')}"
            plot_card(
                L_c,
                B_c,
                Z,
                title,
                cards_dir / f"theta{int(theta)}_{col_key}.png",
            )

    # Sanity vs cone75
    Z75p = maps[(75.0, "phys_invsq")]
    Z75s = maps[(75.0, "sheet_invsq")]
    ref_p = ref["cone_phys_invsq"]
    ref_s = ref["cone_sheet_invsq"]
    diff_p = float(np.nanmax(np.abs(Z75p - ref_p)))
    diff_s = float(np.nanmax(np.abs(Z75s - ref_s)))
    print("\n=== θ=75° sanity vs cone75_maps.npz ===", flush=True)
    print(
        f"  phys_invsq:  ours [{Z75p.min():.6f},{Z75p.max():.6f}] "
        f"ref [{ref_p.min():.6f},{ref_p.max():.6f}] max|Δ|={diff_p:.3e}",
        flush=True,
    )
    print(
        f"  sheet_invsq: ours [{Z75s.min():.6f},{Z75s.max():.6f}] "
        f"ref [{ref_s.min():.6f},{ref_s.max():.6f}] max|Δ|={diff_s:.3e}",
        flush=True,
    )

    # ---- Primary 4×4 grid ----
    style()
    # ~2500 px wide @ 120 dpi → ~20.8 in; height ~ 4*2.6 + margins ≈ 12
    fig_w, fig_h = 22.0, 14.5
    fig, axes = plt.subplots(
        4,
        4,
        figsize=(fig_w, fig_h),
        subplot_kw={"projection": "mollweide"},
        dpi=120,
    )
    pcm = None
    for i, theta in enumerate(THETAS):
        for j, (col_key, _, col_title) in enumerate(COL_KEYS):
            ax = axes[i, j]
            title = col_title if i == 0 else None
            pcm = plot_panel(
                ax,
                L_c,
                B_c,
                maps[(theta, col_key)],
                title=title,
                xticklabels=(i == 3),
            )
        # row label on left of first column via fig text
    # Row labels in figure coords
    for i, theta in enumerate(THETAS):
        # get axes position
        pos = axes[i, 0].get_position()
        y = (pos.y0 + pos.y1) / 2
        fig.text(
            0.012,
            y,
            rf"$\theta={int(theta)}^\circ$",
            rotation=90,
            va="center",
            ha="center",
            fontsize=13,
            fontweight="bold",
        )

    fig.suptitle(
        r"Cone-averaged $A/A_i$ Mollweides at Migkas cone sizes"
        f"  |  magma [{VMIN:.2f}, {VMAX:.2f}]  |  ℓ increases LEFT",
        fontsize=14,
        y=0.98,
    )
    fig.subplots_adjust(
        left=0.05, right=0.90, bottom=0.04, top=0.93, wspace=0.06, hspace=0.18
    )
    cax = fig.add_axes([0.915, 0.12, 0.015, 0.72])
    cbar = fig.colorbar(pcm, cax=cax)
    cbar.set_label(r"$A/A_i$", fontsize=12)
    ticks = np.linspace(VMIN, VMAX, 8)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t:.2f}" for t in ticks])

    grid4 = OUT / "GRID_4x4_cone_thetas.png"
    fig.savefig(grid4, dpi=120, facecolor="white")
    plt.close(fig)
    print(f"\n4×4 → {grid4}", flush=True)

    # ---- Bonus 2×2 sheet invsq ----
    fig2, axes2 = plt.subplots(
        2,
        2,
        figsize=(14.0, 9.5),
        subplot_kw={"projection": "mollweide"},
        dpi=120,
    )
    for idx, theta in enumerate(THETAS):
        r, c = divmod(idx, 2)
        ax = axes2[r, c]
        pcm2 = plot_panel(
            ax,
            L_c,
            B_c,
            maps[(theta, "sheet_invsq")],
            title=rf"$\theta={int(theta)}^\circ$ — Sheet-rotated invsq",
            xticklabels=True,
        )
    fig2.suptitle(
        r"Ours like Migkas Fig. 8 — sheet-rotated $A/A_i$ (invsq) cone averages"
        f"  |  magma [{VMIN:.2f}, {VMAX:.2f}]",
        fontsize=13,
        y=0.98,
    )
    fig2.subplots_adjust(
        left=0.04, right=0.88, bottom=0.05, top=0.90, wspace=0.08, hspace=0.22
    )
    cax2 = fig2.add_axes([0.90, 0.15, 0.02, 0.65])
    cbar2 = fig2.colorbar(pcm2, cax=cax2)
    cbar2.set_label(r"$A/A_i$", fontsize=11)
    cbar2.set_ticks(ticks)
    grid2 = OUT / "GRID_ours_like_fig8_sheet_invsq.png"
    fig2.savefig(grid2, dpi=120, facecolor="white")
    plt.close(fig2)
    print(f"2×2 → {grid2}", flush=True)

    # ---- npz ----
    npz_data = {
        "L_centers": L_c,
        "B_centers": B_c,
        "thetas_deg": np.array(THETAS),
        "vmin": VMIN,
        "vmax": VMAX,
    }
    for theta in THETAS:
        for col_key, _, _ in COL_KEYS:
            npz_data[f"cone_th{int(theta)}_{col_key}"] = maps[(theta, col_key)]
    npz_path = OUT / "cone_thetas_maps.npz"
    np.savez_compressed(npz_path, **npz_data)
    print(f"npz → {npz_path}", flush=True)

    # ---- NOTE ----
    lines = [
        "# Cone-averaged A/A_i at Migkas cone sizes (θ = 90°, 75°, 60°, 45°)",
        "",
        f"**Created:** {now_local()}",
        "",
        "## Deliverables",
        f"- `{grid4.name}` — 4×4 grid (rows θ large→small; cols Physics invsq | Sheet invsq | Physics hratio | Sheet hratio)",
        f"- `{grid2.name}` — 2×2 Migkas-style (sheet-rotated invsq only)",
        f"- `{npz_path.name}` — 16 cone maps + L/B centers",
        "- `cards/` — optional individual Mollweide cards",
        "",
        "## Method",
        "- Full-sky maps from `fullsky_AB/fullsky_AB_maps.npz`",
        "- Cone geometry: same 72×35 L/B centers as `cone75_maps.npz`",
        "- For each θ: pixels with n·n_center ≥ cos(θ), weight-average, then mean-renorm → A/A_i",
        "- Plot: Migkas Galactic (ℓ increases LEFT), magma [0.85, 1.20]",
        "",
        "## θ=75° sanity vs cone75_maps.npz",
        f"- phys_invsq max|Δ| = {diff_p:.3e}",
        f"- sheet_invsq max|Δ| = {diff_s:.3e}",
        f"- ours phys [{Z75p.min():.6f}, {Z75p.max():.6f}] vs ref [{ref_p.min():.6f}, {ref_p.max():.6f}]",
        f"- ours sheet [{Z75s.min():.6f}, {Z75s.max():.6f}] vs ref [{ref_s.min():.6f}, {ref_s.max():.6f}]",
        "",
        "## Ranges (min, max) for all 16 maps",
        "",
        "| θ | Physics invsq | Sheet invsq | Physics hratio | Sheet hratio |",
        "|---|---------------|-------------|----------------|--------------|",
    ]
    for theta in THETAS:
        cells = []
        for col_key, _, _ in COL_KEYS:
            mn, mx = ranges[(theta, col_key)]
            cells.append(f"[{mn:.3f}, {mx:.3f}]")
        lines.append(f"| {int(theta)}° | " + " | ".join(cells) + " |")
    lines += [
        "",
        "## Morphology vs θ",
        "- **θ=90°** (largest): smoothest, lowest contrast — nearly whole-sky average washes features.",
        "- **θ=75°**: Migkas Fig. 8 matching size; moderate contrast.",
        "- **θ=60°**: sharper hot/cold structure.",
        "- **θ=45°** (smallest): highest contrast / sharpest features (least averaging).",
        "",
        "Same trend as Migkas Fig. 8: smaller cones → more contrast.",
        "",
        f"Done in {time.time()-t0:.1f}s.",
    ]
    note_path = OUT / "NOTE.md"
    note_path.write_text("\n".join(lines) + "\n")
    print(f"NOTE → {note_path}", flush=True)
    print(f"\nDone in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()

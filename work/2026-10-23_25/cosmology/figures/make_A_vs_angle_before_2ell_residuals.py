#!/usr/bin/env python3
"""Make the before-2ell A/A_all-versus-angle plot with residuals."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SOURCE = Path("/workspace/cosmology/map_compare/amp_mean_param_noise/best_fit_2ell_multiplier.npz")
FIGURE_DIR = Path("/workspace/cosmology/figures")
COPY_DIR = Path("/workspace/cosmology/map_compare/amp_mean_param_noise")
FILENAME = "A_vs_angle_before_2ell_with_residuals.png"
SIDECAR_NAME = "A_vs_angle_before_2ell_with_residuals.json"
SCRIPT_PATH = FIGURE_DIR / "make_A_vs_angle_before_2ell_residuals.py"


def periodic_curve(x_deg: np.ndarray, y: np.ndarray, n: int = 1441) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate one periodic cycle so the model is rendered as a smooth curve."""
    order = np.argsort(x_deg)
    x = np.asarray(x_deg)[order]
    values = np.asarray(y)[order]
    x_ext = np.r_[x - 360.0, x, x + 360.0]
    y_ext = np.tile(values, 3)
    grid = np.linspace(0.0, 360.0, n)
    return grid, np.interp(grid, x_ext, y_ext)


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    COPY_DIR.mkdir(parents=True, exist_ok=True)

    with np.load(SOURCE) as data:
        ell = np.asarray(data["ell_m"], dtype=float)
        measured = np.asarray(data["A_m"], dtype=float)
        smooth = np.asarray(data["A_before"], dtype=float)
        residual = np.asarray(data["resid_before"], dtype=float)
        r2 = float(data["R2_before"])
        rmse = float(data["RMSE_before"])
        amp1 = float(data["Amp1"])
        phi1 = float(data["phi1"])

    grid, smooth_grid = periodic_curve(ell, smooth)
    gold = "#c9961a"
    charcoal = "#202124"
    grid_color = "#d9dde3"

    plt.rcParams.update({
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "legend.fontsize": 10,
        "figure.dpi": 120,
        "savefig.dpi": 300,
    })
    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(10.5, 7.6),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0], "hspace": 0.18},
    )

    # Top: measured points and only the requested before-2ell smooth model.
    ax_top.plot(
        grid,
        smooth_grid,
        color=gold,
        lw=3.0,
        solid_capstyle="round",
        label="Amp×sin smooth (before 2ℓ)",
        zorder=2,
    )
    ax_top.scatter(
        ell,
        measured,
        s=30,
        color="black",
        edgecolors="white",
        linewidths=0.45,
        label="Migkas data",
        zorder=3,
    )
    ax_top.set_ylabel("A/A_all")
    ax_top.set_title("A/A_all vs angle — before rotation / 2ℓ", loc="left", pad=10, fontweight="bold")
    ax_top.text(
        0.99,
        1.015,
        "Amp×sin smooth model; before 2ℓ correction",
        transform=ax_top.transAxes,
        ha="right",
        va="bottom",
        color="#555b66",
        fontsize=10,
    )
    ax_top.legend(
        loc="upper left",
        frameon=True,
        framealpha=0.94,
        facecolor="white",
        edgecolor="#c7ccd4",
        title=(f"R²_before = {r2:.4f}   RMSE_before = {rmse:.4f}\n"
               f"Amp1 = {amp1:.4f}   phi1 = {phi1:.4f} rad"),
        title_fontsize=9.5,
    )

    # Bottom: supplied before-2ell residuals and zero reference.
    ax_bottom.axhline(0.0, color=charcoal, lw=1.0, ls="--", alpha=0.8, zorder=1)
    ax_bottom.plot(ell, residual, color="#3d566e", lw=1.2, alpha=0.7, zorder=2)
    ax_bottom.scatter(
        ell,
        residual,
        s=24,
        color="#3d566e",
        edgecolors="white",
        linewidths=0.4,
        zorder=3,
    )
    ax_bottom.set_ylabel("Residual")
    ax_bottom.set_xlabel("Galactic longitude / scan angle (deg)")
    ax_bottom.set_title(f"Residuals (A_m − A_before) — before 2ℓ; RMSE_before = {rmse:.4f}", loc="left", pad=8)

    # Shared clean axes and requested 45-degree ticks.
    ticks = np.arange(0.0, 361.0, 45.0)
    for ax in (ax_top, ax_bottom):
        ax.set_xlim(0.0, 360.0)
        ax.set_xticks(ticks)
        ax.grid(axis="y", color=grid_color, lw=0.8, alpha=0.7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    ax_bottom.set_xticklabels([f"{int(t)}°" for t in ticks])

    # Use compact, symmetric residual limits with a little breathing room.
    resid_limit = max(0.1, 1.08 * float(np.nanmax(np.abs(residual))))
    ax_bottom.set_ylim(-resid_limit, resid_limit)
    top_min = min(float(np.nanmin(measured)), float(np.nanmin(smooth_grid)))
    top_max = max(float(np.nanmax(measured)), float(np.nanmax(smooth_grid)))
    top_pad = max(0.015, 0.08 * (top_max - top_min))
    ax_top.set_ylim(top_min - top_pad, top_max + top_pad)

    fig.suptitle(
        "Migkas A/A_all compared with the smooth model (before rotation / 2ℓ)",
        y=0.995,
        fontsize=14,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.91, bottom=0.10, left=0.10, right=0.98, hspace=0.18)

    primary = FIGURE_DIR / FILENAME
    copied = COPY_DIR / FILENAME
    fig.savefig(primary, bbox_inches="tight", facecolor="white")
    fig.savefig(copied, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    sidecar = {
        "R2_before": r2,
        "RMSE_before": rmse,
        "Amp1": amp1,
        "phi1": phi1,
        "source_npz": str(SOURCE),
        "figure_path": str(primary),
        "figure_copy_path": str(copied),
        "script_path": str(SCRIPT_PATH),
    }
    (FIGURE_DIR / SIDECAR_NAME).write_text(json.dumps(sidecar, indent=2) + "\n")
    (COPY_DIR / SIDECAR_NAME).write_text(json.dumps({**sidecar, "sidecar_path": str(COPY_DIR / SIDECAR_NAME)}, indent=2) + "\n")

    print(json.dumps({**sidecar, "sidecar_path": str(FIGURE_DIR / SIDECAR_NAME)}, indent=2))


if __name__ == "__main__":
    main()

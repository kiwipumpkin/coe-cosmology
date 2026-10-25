#!/usr/bin/env python3
"""Make Chris-format axisymmetric before/after rotation plot with residuals."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

CSV_AF = Path("/workspace/cosmology/map_compare/amp_mean_param_noise/A_vs_angle_ideal_plus_chris_rotation.csv")
CSV_MIGKAS = Path("/workspace/cosmology/migkas/fig6_digitized.csv")
FIGURE_DIR = Path("/workspace/cosmology/figures")
COPY_DIR = Path("/workspace/cosmology/map_compare/amp_mean_param_noise")
FILENAME = "A_vs_angle_axisymmetric_before_rotation_with_residuals.png"
JSON_NAME = "A_vs_angle_axisymmetric_before_rotation_with_residuals.json"
SCRIPT_PATH = FIGURE_DIR / "make_A_vs_angle_axisymmetric_before_rotation.py"

BLUE = "#1f77b4"
GOLD = "#c49a17"
BLACK = "#111111"
GRAY = "#777777"
RED = "#d62728"


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
            "legend.fontsize": 9,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 1.1,
            "mathtext.fontset": "dejavusans",
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )


def read_migkas(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {
        "ell": np.array([float(row["l_deg"]) for row in rows]),
        "A": np.array([float(row["A_over_Aall"]) for row in rows]),
        "err_low": np.array([float(row["A_err_low"]) for row in rows]),
        "err_high": np.array([float(row["A_err_high"]) for row in rows]),
    }


def periodic_interp(x_deg: np.ndarray, y: np.ndarray, query_deg: np.ndarray) -> np.ndarray:
    """Interpolate a 0..360 curve with periodic wrapping at the boundaries."""
    order = np.argsort(x_deg)
    x = np.asarray(x_deg, dtype=float)[order]
    values = np.asarray(y, dtype=float)[order]
    x_ext = np.concatenate((x - 360.0, x, x + 360.0))
    y_ext = np.tile(values, 3)
    return np.interp(np.mod(np.asarray(query_deg, dtype=float), 360.0), x_ext, y_ext)


def r2_rmse(observed: np.ndarray, modeled: np.ndarray) -> tuple[float, float]:
    residual = np.asarray(observed) - np.asarray(modeled)
    r2 = 1.0 - float(np.sum(residual**2)) / float(np.sum((observed - np.mean(observed))**2))
    rmse = float(np.sqrt(np.mean(residual**2)))
    return r2, rmse


def main() -> None:
    style()
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    COPY_DIR.mkdir(parents=True, exist_ok=True)

    # AF_smooth is the axisymmetric before-rotation curve; A_rotated_smooth is after.
    # No Amp×sin or 2ell model is used here.
    af_data = np.genfromtxt(CSV_AF, delimiter=",", names=True)
    ell_af = np.asarray(af_data["ell_deg"], dtype=float)
    AF_smooth = np.asarray(af_data["AF_smooth"], dtype=float)
    A_rotated_smooth = np.asarray(af_data["A_rotated_smooth"], dtype=float)
    mig = read_migkas(CSV_MIGKAS)
    ell = mig["ell"]
    A = mig["A"]
    A_err_low = mig["err_low"]
    A_err_high = mig["err_high"]

    af_before = periodic_interp(ell_af, AF_smooth, ell)
    af_after = periodic_interp(ell_af, A_rotated_smooth, ell)
    residual_before = A - af_before
    correction = af_after - af_before
    r2_before, rmse_before = r2_rmse(A, af_before)
    r2_after, rmse_after = r2_rmse(A, af_after)

    ell_fine = np.linspace(0.0, 360.0, 1441)
    before_fine = periodic_interp(ell_af, AF_smooth, ell_fine)
    after_fine = periodic_interp(ell_af, A_rotated_smooth, ell_fine)
    correction_fine = after_fine - before_fine

    fig = plt.figure(figsize=(10.5, 8.2))
    gs = GridSpec(2, 1, height_ratios=[1, 1], hspace=0.08, figure=fig)
    ax_top = fig.add_subplot(gs[0])
    ax_bottom = fig.add_subplot(gs[1], sharex=ax_top)

    ax_top.plot(
        ell_fine,
        before_fine,
        color=BLUE,
        lw=2.2,
        ls="--",
        label=f"before rotation: axisymmetric AF ($R^2={r2_before:.3f}$)",
        zorder=2,
    )
    ax_top.plot(
        ell_fine,
        after_fine,
        color=GOLD,
        lw=2.8,
        label=f"after sine rotation ($R^2={r2_after:.3f}$)",
        zorder=2,
    )
    ax_top.errorbar(
        ell,
        A,
        yerr=np.vstack((A_err_low, A_err_high)),
        fmt="o",
        ms=4.0,
        color=BLACK,
        ecolor=GRAY,
        elinewidth=0.9,
        capsize=2.0,
        label="Migkas",
        zorder=3,
    )
    ax_top.axhline(1.0, color="#aaaaaa", lw=0.8, ls=":", zorder=1)
    ax_top.set_ylabel(r"$A/A_{\mathrm{all}}$")
    ax_top.set_ylim(0.75, 1.25)
    ax_top.set_title("Axisymmetric AF + sine rotation — A vs Galactic longitude", pad=8)
    ax_top.legend(loc="upper right", framealpha=0.92)

    ax_bottom.plot(
        ell_fine,
        correction_fine,
        color=RED,
        lw=2.8,
        label=r"rotation correction: $A_{\rm after}-A_{\rm before}$",
        zorder=2,
    )
    ax_bottom.errorbar(
        ell,
        residual_before,
        yerr=np.vstack((A_err_low, A_err_high)),
        fmt="o",
        ms=4.0,
        color="#666666",
        ecolor="#888888",
        elinewidth=0.9,
        capsize=2.0,
        label="residual before rotation",
        zorder=3,
    )
    ax_bottom.axhline(0.0, color="black", lw=0.9, ls="-", zorder=1)
    ax_bottom.set_ylabel("Migkas − model")
    ax_bottom.set_ylim(-0.20, 0.20)
    ax_bottom.set_title(
        f"Residuals before rotation (RMSE {rmse_before:.3f}; after rotation would be {rmse_after:.3f})",
        pad=8,
    )
    ax_bottom.set_xlabel(r"Galactic longitude $\ell$ [deg]")
    ax_bottom.legend(loc="upper right", framealpha=0.92)

    ticks = np.arange(0.0, 351.0, 50.0)
    for ax in (ax_top, ax_bottom):
        ax.set_xlim(0.0, 360.0)
        ax.set_xticks(ticks)
        ax.grid(True, color="#dddddd", lw=0.8, alpha=0.85)
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.spines["top"].set_visible(True)
        ax.spines["right"].set_visible(True)
    plt.setp(ax_top.get_xticklabels(), visible=False)

    primary = FIGURE_DIR / FILENAME
    copied = COPY_DIR / FILENAME
    fig.savefig(primary, bbox_inches="tight", facecolor="white")
    fig.savefig(copied, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    sidecar = {
        "R2_before": r2_before,
        "RMSE_before": rmse_before,
        "R2_after": r2_after,
        "RMSE_after": rmse_after,
        "note": "Before curve is AF_smooth and after curve is A_rotated_smooth from A_vs_angle_ideal_plus_chris_rotation.csv; bottom red curve is A_rotated_smooth - AF_smooth. This is the axisymmetric rotation analogue, not Amp×sin or 2ell.",
    }
    (FIGURE_DIR / JSON_NAME).write_text(json.dumps(sidecar, indent=2) + "\n")
    (COPY_DIR / JSON_NAME).write_text(json.dumps(sidecar, indent=2) + "\n")
    print(json.dumps({**sidecar, "figure": str(primary), "copy": str(copied), "script": str(SCRIPT_PATH)}, indent=2))


if __name__ == "__main__":
    main()

"""Plot photon energy change needed to cancel the mean-subtracted Doppler residual."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "dtt.json"), encoding="utf-8") as f:
    data = json.load(f)
rows = data["rows"]

OUT = os.path.join(HERE, "kill_energy_pct_0_360.png")
LIMIT_FRAC = 1.0e-5
LIMIT_PCT = 100.0 * LIMIT_FRAC  # 0.001 percent
RED = "#c43d2f"
GREEN = "#238b45"
BLACK = "#202020"

plt.rcParams.update({
    "font.size": 15,
    "axes.titlesize": 17,
    "axes.labelsize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
})


def full_circle(values):
    values = np.asarray(values, dtype=float)
    return np.concatenate((values, values[-2:0:-1], values[:1]))


def pct_fmt(value, _pos):
    if value == 0:
        return "0"
    if abs(value) >= 1:
        return f"{value:g}"
    return f"{value:.3g}"


fig, axes = plt.subplots(3, 1, figsize=(22, 21), sharex=True)
label_positions = {
    "14.7": (274, -30.0),
    "22.28": (202, 205.0),
    "49": (245, -23.5),
}

for ax, row in zip(axes, rows):
    theta = np.asarray(row["theta_toward_deg"], dtype=float)
    phi = np.concatenate((theta, 360.0 - theta[-2:0:-1], np.array([360.0])))
    doppler = np.asarray(row["dtt_doppler"], dtype=float)
    kill_pct = -100.0 * (doppler - doppler.mean())
    y = full_circle(kill_pct)
    dlab = row["D_label"]
    peak = float(np.max(np.abs(y)))
    ymin = float(np.min(y))
    ymax = float(np.max(y))

    # Symmetric-log keeps the +/-0.001 percent CMB band visible.
    ax.set_yscale("symlog", linthresh=0.01, linscale=1.0, base=10)
    ax.axhspan(-LIMIT_PCT, LIMIT_PCT, color=GREEN, alpha=0.18, zorder=0)
    ax.axhline(+LIMIT_PCT, color=GREEN, lw=1.5, ls="--", zorder=1)
    ax.axhline(-LIMIT_PCT, color=GREEN, lw=1.5, ls="--", zorder=1)
    ax.axhline(0, color=BLACK, lw=0.8, alpha=0.75, zorder=1)
    ax.plot(phi, y, color=RED, lw=4.8, solid_capstyle="round", zorder=4)

    pad = max(0.08 * max(abs(ymin), abs(ymax)), 0.01)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 30))
    ax.yaxis.set_major_formatter(FuncFormatter(pct_fmt))
    ax.grid(alpha=0.24, linewidth=0.8, which="both")
    ax.set_ylabel("photon dE/E (%)")
    ax.set_title(
        f"D = {dlab} Gpc    |    peak required |dE/E| = {peak:.3f}%    "
        f"(range {ymin:+.3f}% to {ymax:+.3f}%)",
        loc="left", pad=11, fontweight="bold",
    )

    lx, ly = label_positions[dlab]
    ax.text(
        lx, ly, "required correction\n+ gain / - loss", color=RED,
        fontsize=13.5, ha="left", va="center", zorder=8,
        bbox=dict(boxstyle="round,pad=0.32", fc="white", ec="none", alpha=0.90),
    )
    ax.text(
        0.995, 0.505, "observed limit\n+/- 0.001%", transform=ax.transAxes,
        color=GREEN, fontsize=11.5, ha="right", va="center", zorder=8,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.86),
    )

axes[-1].set_xlabel(
    "great-circle angle phi (degrees) about the CoE axis    "
    "0 deg = toward CoE    |    180 deg = anti-CoE    |    360 deg = toward CoE"
)
fig.suptitle(
    "Photon energy correction required to cancel the kinematic Doppler residual\n"
    "full 0-360 deg great circle; axisymmetric reflection of the theta = 0-180 deg calculation",
    fontsize=21, fontweight="bold", y=0.985,
)
fig.text(
    0.5, 0.945,
    "Thick red: 100 * [-dtt_doppler + mean(dtt_doppler)] = required dE/E (%)    "
    "green band/lines: |dT/T| = 1e-5 = +/- 0.001%",
    ha="center", va="center", fontsize=14,
)
fig.text(
    0.5, 0.018,
    "Positive dE/E means the photon must gain energy; negative means it must lose energy. "
    "Y-axis is symmetric-log so the CMB-limit band is visible.",
    ha="center", va="center", fontsize=12.5, color="#444444",
)
fig.tight_layout(rect=(0.035, 0.045, 0.995, 0.925), h_pad=2.2)
fig.savefig(OUT, dpi=200, facecolor="white")
plt.close(fig)

print("wrote", OUT)
for row in rows:
    dop = np.asarray(row["dtt_doppler"], dtype=float)
    pct = -100.0 * (dop - dop.mean())
    print(
        f"D={row['D_label']} Gpc: peak_abs={np.max(np.abs(pct)):.9f}% "
        f"range={np.min(pct):+.9f}%..{np.max(pct):+.9f}%"
    )

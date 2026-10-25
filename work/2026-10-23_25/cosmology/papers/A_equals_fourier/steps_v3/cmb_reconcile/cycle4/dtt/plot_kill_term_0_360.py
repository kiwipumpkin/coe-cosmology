"""Plot the required CMB kill term over a full great circle (0--360 deg).

The source calculation is tabulated for theta=0..180 deg from the CoE axis.
For an axisymmetric field, phi and 360-phi are identical, so the half-circle
is mirrored about the CoE/anti-CoE line.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "dtt.json")) as f:
    data = json.load(f)
rows = data["rows"]

KILL = "#c43d2f"
DOPPLER = "#1769aa"

plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 15,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
})


def demean(values):
    values = np.asarray(values, dtype=float)
    return values - values.mean()


def full_circle(values):
    """Mirror values tabulated at theta=0..180 onto phi=0..360."""
    values = np.asarray(values, dtype=float)
    # 0..180, then 181..359 from theta=179..1, and the closing 360 point.
    return np.concatenate((values, values[-2:0:-1], values[:1]))


def label(ax, x, y, text, color, *, ha="left", size=12.5):
    ax.text(
        x, y, text, color=color, fontsize=size, ha=ha, va="center", zorder=8,
        linespacing=1.18,
        bbox=dict(boxstyle="round,pad=0.30", fc="white", ec="none", alpha=0.88),
    )

fig, axes = plt.subplots(3, 1, figsize=(18, 19), sharex=True)

# Positions are intentionally in open portions of each panel; no leader lines.
positions = {
    "14.7": {"kill": (218, 0.36), "dop": (218, -0.34)},
    "22.28": {"kill": (215, 2.05), "dop": (215, -1.95)},
    "49": {"kill": (218, 0.43), "dop": (218, -0.43)},
}

for ax, row in zip(axes, rows):
    theta = np.asarray(row["theta_toward_deg"], dtype=float)
    phi = np.concatenate((theta, 360.0 - theta[-2:0:-1], np.array([360.0])))
    doppler = demean(row["dtt_doppler"])
    kill = -doppler
    doppler_full = full_circle(doppler)
    kill_full = full_circle(kill)
    dlab = row["D_label"]
    kill_a2 = -float(row["doppler_a2"])

    ax.plot(phi, kill_full, color=KILL, lw=4.8, solid_capstyle="round",
            zorder=4, label="required kill term")
    ax.plot(phi, doppler_full, color=DOPPLER, lw=1.45, alpha=0.95,
            zorder=3, label="Doppler residual")
    ax.axhline(0, color="black", lw=0.75, alpha=0.72, zorder=1)
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 30))
    ax.grid(alpha=0.24, linewidth=0.8)
    ax.set_ylabel(r"$\Delta T/T$  (mean-subtracted)")

    if dlab == "22.28":
        ylim = (-3.05, 3.05)
    elif dlab == "14.7":
        ylim = (-0.52, 0.52)
    else:
        ylim = (-0.62, 0.62)
    ax.set_ylim(*ylim)

    ax.set_title(
        "D = %s Gpc     required kill $a_2$ = %+.3f" % (dlab, kill_a2),
        loc="left", pad=10, fontweight="bold",
    )

    xk, yk = positions[dlab]["kill"]
    xd, yd = positions[dlab]["dop"]
    label(ax, xk, yk, "required kill term\n$-\\,\\Delta T/T_{\\rm Doppler}$", KILL)
    label(ax, xd, yd, "Doppler residual\n(thin)", DOPPLER)

axes[-1].set_xlabel(
    r"great-circle angle $\phi$ about the CoE axis (degrees)    "
    r"$0^\circ$ = toward CoE    ·    $180^\circ$ = anti-CoE    ·    "
    r"$360^\circ$ = toward CoE"
)

fig.suptitle(
    "Required CMB kill term vs angle: full 0–360° great circle\n"
    "Chris's CoE model · axisymmetric mirror of the θ = 0–180° calculation",
    fontsize=20, fontweight="bold", y=0.985,
)
fig.text(
    0.5, 0.948,
    r"Thick red: required $-\,\Delta T/T_{\rm Doppler}$ (mean-subtracted)   ·   "
    r"thin blue: Doppler residual",
    ha="center", va="center", fontsize=14,
)
fig.tight_layout(rect=(0.03, 0.035, 0.995, 0.925), h_pad=2.0)

out = os.path.join(HERE, "kill_term_0_360.png")
fig.savefig(out, dpi=180, facecolor="white")
print("wrote", out)
for row in rows:
    print("D=%s Gpc: kill a2=%+.9f" %
          (row["D_label"], -float(row["doppler_a2"])))

"""Plot the CMB term required to cancel the kinematic Doppler residual.

The plotted kill term is -demean(Doppler), so Doppler + kill is identically
zero (up to floating-point roundoff).  All angles use theta=0 toward CoE.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(HERE, "dtt.json")))
ROWS = DATA["rows"]

KILL = "#c43d2f"       # emphasized required term
DOPPLER = "#1769aa"    # thin kinematic residual
RESID = "#666666"

plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 15,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
})


def demean(y):
    y = np.asarray(y, dtype=float)
    return y - y.mean()


def boxed_text(ax, x, y, text, color, *, ha="left", va="center", size=12.5):
    """Place an in-panel label without a leader/pointer line."""
    ax.text(
        x, y, text, color=color, fontsize=size, ha=ha, va=va,
        zorder=8, linespacing=1.18,
        bbox=dict(boxstyle="round,pad=0.30", fc="white", ec="none", alpha=0.88),
    )

fig, axes = plt.subplots(3, 1, figsize=(16, 18.5), sharex=True)

# Hand-picked open-space label locations; labels deliberately have no leaders.
label_positions = {
    "14.7": {
        "kill": (112, 0.34),
        "dop": (42, -0.31),
        "resid": (146, 0.035),
    },
    "22.28": {
        "kill": (14, 2.25),
        "dop": (14, -2.10),
        "resid": (126, 0.075),
    },
    "49": {
        "kill": (108, 0.40),
        "dop": (108, 0.24),
        "resid": (145, 0.045),
    },
}

for ax, row in zip(axes, ROWS):
    th = np.asarray(row["theta_toward_deg"], dtype=float)
    doppler = demean(row["dtt_doppler"])
    kill = -doppler
    residual = doppler + kill
    dlab = row["D_label"]
    kill_a2 = -float(row["doppler_a2"])

    lim = max(float(np.max(np.abs(doppler))), float(np.max(np.abs(kill))))
    # Extra headroom makes the direct labels sit in open space.
    if dlab == "22.28":
        ylim = (-3.05, 3.05)
    elif dlab == "14.7":
        ylim = (-0.52, 0.52)
    else:
        ylim = (-0.62, 0.62)

    ax.plot(th, kill, color=KILL, lw=4.6, solid_capstyle="round", zorder=4)
    ax.plot(th, doppler, color=DOPPLER, lw=1.45, alpha=0.95, zorder=3)
    ax.plot(th, residual, color=RESID, lw=1.8, ls=(0, (2, 2)), zorder=2)
    ax.axhline(0, color="black", lw=0.75, alpha=0.7, zorder=1)
    ax.set_xlim(0, 180)
    ax.set_ylim(*ylim)
    ax.set_xticks([0, 30, 60, 90, 120, 150, 180])
    ax.grid(alpha=0.24, linewidth=0.8)
    ax.set_ylabel(r"$\Delta T/T$  (mean-subtracted)")

    ax.set_title(
        "D = %s Gpc     required kill $a_2$ = %+.3f" % (dlab, kill_a2),
        loc="left", pad=10, fontweight="bold",
    )

    # Direct labels in open space, with no leader/pointer lines.
    xk, yk = label_positions[dlab]["kill"]
    xd, yd = label_positions[dlab]["dop"]
    xr, yr = label_positions[dlab]["resid"]
    boxed_text(ax, xk, yk, "required kill term\n$-\\,\\Delta T/T_{\\rm Doppler}$", KILL)
    boxed_text(ax, xd, yd, "Doppler residual\n(thin)", DOPPLER)
    boxed_text(ax, xr, yr, "after cancel\n$\\simeq 0$", RESID, ha="center", size=11.5)

    if dlab == "22.28":
        boxed_text(ax, 52, -2.72, "turnaround spike is canceled point-by-point", "0.30", size=11.5)

axes[-1].set_xlabel(
    r"angle $\theta$ from CoE axis (degrees)    "
    r"$0^\circ$ = toward CoE    ·    $180^\circ$ = away from CoE"
)

fig.suptitle(
    "Required CMB kill term vs angle\n"
    "Chris's CoE model: cancel the kinematic Doppler residual on the CMB emission sphere",
    fontsize=19, fontweight="bold", y=0.985,
)
fig.text(
    0.5, 0.948,
    r"Thick red: required $-\,\Delta T/T_{\rm Doppler}$ (mean-subtracted)   ·   "
    r"thin blue: Doppler residual   ·   gray dotted: Doppler + kill",
    ha="center", va="center", fontsize=13.5,
)
fig.tight_layout(rect=(0.03, 0.035, 0.995, 0.925), h_pad=2.0)

out = os.path.join(HERE, "kill_term_vs_theta.png")
fig.savefig(out, dpi=180, facecolor="white")
print("wrote", out)
for row in ROWS:
    print("D=%s Gpc: kill a2=%+.9f" % (row["D_label"], -float(row["doppler_a2"])))

"""0–360° kill-term plots for a few recombination-track points (incl. headline-like)."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(HERE, "track.json")))
rows_all = [r for r in data["rows"] if "kill_a2" in r]

# Pick track points: near shell-min, mid, static-ish, headline
WANT = [4.5, 14.7, 22.28, 49.17]
rows = []
for w in WANT:
    best = min(rows_all, key=lambda r: abs(r["D"] - w))
    rows.append(best)

KILL = "#c43d2f"
DOPPLER = "#1769aa"

plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})


def demean(v):
    v = np.asarray(v, dtype=float)
    return v - v.mean()


def full_circle(values):
    values = np.asarray(values, dtype=float)
    return np.concatenate((values, values[-2:0:-1], values[:1]))


fig, axes = plt.subplots(len(rows), 1, figsize=(18, 5.2 * len(rows)), sharex=True)
if len(rows) == 1:
    axes = [axes]

# label positions in open panel regions (no leaders)
positions = {
    4.5: {"kill": (210, 0.085), "dop": (210, -0.085)},
    14.7: {"kill": (210, 0.13), "dop": (210, -0.13)},
    22.28: {"kill": (210, 0.13), "dop": (210, -0.13)},
    49.17: {"kill": (210, 0.12), "dop": (210, -0.12)},
}

for ax, row in zip(axes, rows):
    theta = np.asarray(row["theta_toward_deg"], dtype=float)
    phi = np.concatenate((theta, 360.0 - theta[-2:0:-1], np.array([360.0])))
    doppler = demean(row["dtt_doppler"])
    kill = -doppler
    doppler_full = full_circle(doppler)
    kill_full = full_circle(kill)
    D = row["D"]
    kill_a2 = float(row["kill_a2"])
    R_see = float(row["R_see"])
    D_rec = float(row["D_rec"])

    ax.plot(phi, kill_full, color=KILL, lw=4.6, solid_capstyle="round", zorder=4)
    ax.plot(phi, doppler_full, color=DOPPLER, lw=1.4, alpha=0.95, zorder=3)
    ax.axhline(0, color="black", lw=0.75, alpha=0.7, zorder=1)
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 30))
    ax.grid(alpha=0.24, linewidth=0.8)
    ax.set_ylabel(r"$\Delta T/T$  (mean-subtracted)")

    amp = float(np.max(np.abs(kill_full)))
    pad = 0.18 * amp + 0.01
    ax.set_ylim(-amp - pad, amp + pad)

    tag = "headline-like" if abs(D - 49.17) < 0.2 else ("near shell-min" if D < 5 else "")
    title = (
        r"$D_{{\rm now}}$ = {:.2f} Gpc   ·   $R_{{\rm see}}$ = {:.3f} Gpc   ·   "
        r"$D_{{\rm rec}}$ = {:.4f} Gpc   ·   kill $a_2$ = {:+.4f}"
    ).format(D, R_see, D_rec, kill_a2)
    if tag:
        title += "   ·   " + tag
    ax.set_title(title, loc="left", pad=10, fontweight="bold")

    # nearest position key
    key = min(positions.keys(), key=lambda k: abs(k - D))
    xk, yk = positions[key]["kill"]
    xd, yd = positions[key]["dop"]
    # scale label y to panel
    yk = np.sign(yk) * 0.55 * amp if amp > 0 else yk
    yd = np.sign(yd) * 0.55 * amp if amp > 0 else yd
    ax.text(xk, yk, "required kill term\n$-\\,\\Delta T/T_{\\rm Doppler}$", color=KILL,
            fontsize=12.5, ha="left", va="center", zorder=8,
            bbox=dict(boxstyle="round,pad=0.30", fc="white", ec="none", alpha=0.88))
    ax.text(xd, yd, "Doppler residual\n(thin)", color=DOPPLER,
            fontsize=12.5, ha="left", va="center", zorder=8,
            bbox=dict(boxstyle="round,pad=0.30", fc="white", ec="none", alpha=0.88))

axes[-1].set_xlabel(
    r"great-circle angle $\phi$ about the CoE axis (degrees)    "
    r"$0^\circ$ = toward CoE    ·    $180^\circ$ = anti-CoE    ·    $360^\circ$ = toward CoE"
)
fig.suptitle(
    "Required CMB kill term on the observer-centred $R_{\\rm see}$ sphere\n"
    "Chris recombination track · axisymmetric 0–360° · NOT fixed R = 14.3 Gpc",
    fontsize=19, fontweight="bold", y=0.995,
)
fig.text(
    0.5, 0.955,
    r"Thick red: required $-\,\Delta T/T_{\rm Doppler}$ (mean-subtracted)   ·   thin blue: Doppler residual",
    ha="center", va="center", fontsize=13.5,
)
fig.tight_layout(rect=(0.03, 0.03, 0.995, 0.93), h_pad=1.8)
out = os.path.join(HERE, "kill_term_0_360_track.png")
fig.savefig(out, dpi=180, facecolor="white")
plt.close(fig)
print("wrote", out)
for row in rows:
    print("D=%.2f R_see=%.3f kill_a2=%+.6f peak_pct=%.3f"
          % (row["D"], row["R_see"], row["kill_a2"], row["kill_energy_pct_peak"]))

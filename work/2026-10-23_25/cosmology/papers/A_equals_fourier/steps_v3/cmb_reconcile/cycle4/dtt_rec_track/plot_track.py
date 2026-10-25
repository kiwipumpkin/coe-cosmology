"""Plot kill a2 and peak |dE/E|% along the recombination track."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
data = json.load(open(os.path.join(HERE, "track.json")))
rows = [r for r in data["rows"] if "kill_a2" in r]

PASS = 1e-5
RED = "#c43d2f"
BLUE = "#1769aa"
GREEN = "#238b45"
ORANGE = "#d94801"
BLACK = "#202020"

plt.rcParams.update({
    "font.size": 15,
    "axes.titlesize": 17,
    "axes.labelsize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
})

D = np.array([r["D"] for r in rows])
D_rec = np.array([r["D_rec"] for r in rows])
R_see = np.array([r["R_see"] for r in rows])
kill_a2 = np.array([r["kill_a2"] for r in rows])
peak_pct = np.array([r["kill_energy_pct_peak"] for r in rows])
stat_a2 = np.array([abs(r["static_a2"]) for r in rows])

# ---- Figure 1: kill a2 and peak |dE/E| vs D_rec (twin of D_now) ----
fig, axes = plt.subplots(2, 1, figsize=(16, 14), sharex=False)

ax = axes[0]
ax.plot(D_rec * 1e3, kill_a2, "o-", color=RED, lw=3.2, ms=9, zorder=4, label=r"kill $a_2 = -a_2^{\rm Doppler}$")
ax.plot(D_rec * 1e3, stat_a2, "s--", color=BLUE, lw=2.0, ms=7, zorder=3, label=r"$|a_2^{\rm static}|$")
ax.axhline(PASS, color=GREEN, lw=2.0, ls="--", zorder=2)
ax.set_yscale("log")
ax.set_ylabel(r"$|a_2|$")
ax.set_title("Required kill-term quadrupole along Chris's recombination track", loc="left", fontweight="bold", pad=10)
ax.grid(alpha=0.28, which="both")
# labels in open space, no leaders
ax.text(0.03, 0.22, r"CMB limit $10^{-5}$", transform=ax.transAxes, color=GREEN, fontsize=13,
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.9))
ax.text(0.55, 0.82, r"kill $a_2$ (Doppler)", transform=ax.transAxes, color=RED, fontsize=14,
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.9))
ax.text(0.55, 0.70, r"$|a_2^{\rm static}|$", transform=ax.transAxes, color=BLUE, fontsize=14,
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.9))
ymin = min(kill_a2.min(), stat_a2.min(), PASS) * 0.4
ymax = max(kill_a2.max(), stat_a2.max()) * 3
ax.set_ylim(ymin, ymax)
# secondary top axis: D_now
ax2 = ax.twiny()
ax2.set_xlim(D_rec.min() * 1e3, D_rec.max() * 1e3)
# map D_rec ticks to D_now labels at same positions
xticks = D_rec * 1e3
ax.set_xticks(xticks)
ax.set_xticklabels(["%.2f" % (x) for x in D_rec * 1e3], rotation=40, ha="right", fontsize=11)
ax2.set_xticks(xticks)
ax2.set_xticklabels(["%.1f" % d for d in D], rotation=40, ha="left", fontsize=11)
ax.set_xlabel(r"$D_{\rm rec}$ (Mpc)")
ax2.set_xlabel(r"$D_{\rm now}$ (Gpc)")
ax.text(0.98, 0.08, "no point reaches $10^{-5}$", transform=ax.transAxes, color=BLACK,
        fontsize=13, ha="right",
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.9))

ax = axes[1]
ax.plot(D_rec * 1e3, peak_pct, "o-", color=ORANGE, lw=3.2, ms=9, zorder=4)
ax.axhline(100 * PASS, color=GREEN, lw=2.0, ls="--", zorder=2)  # 0.001%
ax.set_yscale("log")
ax.set_ylabel(r"peak $|{\rm d}E/E|$ required (%)")
ax.set_xlabel(r"$D_{\rm rec}$ (Mpc)")
ax.set_title(r"Peak photon energy correction $100\times|\mathrm{kill}|$ vs $D_{\rm rec}$", loc="left", fontweight="bold", pad=10)
ax.grid(alpha=0.28, which="both")
ax.set_xticks(xticks)
ax.set_xticklabels(["%.2f" % (x) for x in D_rec * 1e3], rotation=40, ha="right", fontsize=11)
ax.text(0.03, 0.18, r"CMB ±0.001%", transform=ax.transAxes, color=GREEN, fontsize=13,
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.9))
ax.text(0.55, 0.82, r"peak |dE/E| on 0–180°", transform=ax.transAxes, color=ORANGE, fontsize=14,
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.9))
# annotate R_see
Rmean = float(np.mean(R_see))
ax.text(0.98, 0.08, r"$R_{\rm see}\approx%.2f$ Gpc (nearly flat on track)" % Rmean,
        transform=ax.transAxes, color=BLACK, fontsize=13, ha="right",
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.9))

fig.suptitle(
    "Recombination track (NOT fixed R=14.3):  $R_{\\rm see}=c\\,\\Delta t$,  "
    "$D_{\\rm rec}/D_{\\rm now}=T_0/3000\\,{\\rm K}$\n"
    "compact core re-solved to cluster $a_2=0.084$  ·  $\\Omega_s=0.049$ shell on $R_{\\rm see}$",
    fontsize=18, fontweight="bold", y=0.995,
)
fig.tight_layout(rect=(0.03, 0.02, 0.99, 0.92), h_pad=2.4)
out1 = os.path.join(HERE, "kill_a2_peak_vs_Drec.png")
fig.savefig(out1, dpi=180, facecolor="white")
plt.close(fig)
print("wrote", out1)

# ---- Figure 2: same vs R_see (shows flat R) and vs D_now ----
fig, axes = plt.subplots(1, 2, figsize=(18, 8))
ax = axes[0]
ax.plot(D, kill_a2, "o-", color=RED, lw=3.2, ms=9)
ax.axhline(PASS, color=GREEN, lw=2.0, ls="--")
ax.set_yscale("log")
ax.set_xlabel(r"$D_{\rm now}$ (Gpc)")
ax.set_ylabel(r"kill $a_2$")
ax.set_title("kill $a_2$ vs $D_{\\rm now}$", loc="left", fontweight="bold")
ax.grid(alpha=0.28, which="both")
ax.text(0.05, 0.2, r"CMB $10^{-5}$", transform=ax.transAxes, color=GREEN, fontsize=13,
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.9))
ax.text(0.55, 0.85, "headline D≈49", transform=ax.transAxes, color=RED, fontsize=13,
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.9))

ax = axes[1]
ax.plot(R_see, kill_a2, "o-", color=RED, lw=3.2, ms=9)
ax.axhline(PASS, color=GREEN, lw=2.0, ls="--")
ax.set_yscale("log")
ax.set_xlabel(r"$R_{\rm see}$ (Gpc)")
ax.set_ylabel(r"kill $a_2$")
ax.set_title(r"kill $a_2$ vs $R_{\rm see}$ (track is nearly vertical)", loc="left", fontweight="bold")
ax.grid(alpha=0.28, which="both")
ax.text(0.05, 0.2, r"$R_{\rm see}\approx c\,t_{\rm now}(1-T_0/3000)$", transform=ax.transAxes,
        color=BLACK, fontsize=12,
        bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="none", alpha=0.9))
span = max(R_see.max() - R_see.min(), 0.02)
ax.set_xlim(R_see.min() - 0.3 * span - 0.01, R_see.max() + 0.3 * span + 0.01)

fig.suptitle("Track diagnostics: Doppler kill quadrupole vs today's CoE distance and vs $R_{\\rm see}$",
             fontsize=17, fontweight="bold", y=0.98)
fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.92))
out2 = os.path.join(HERE, "kill_a2_vs_Dnow_Rsee.png")
fig.savefig(out2, dpi=180, facecolor="white")
plt.close(fig)
print("wrote", out2)

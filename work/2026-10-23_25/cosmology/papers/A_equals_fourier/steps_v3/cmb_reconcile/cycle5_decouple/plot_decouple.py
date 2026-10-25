"""Full-size plots: |ΔH/H| vs |ΔT/T| components from the same Φ. Labels in open space, no leader lines."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(HERE, "decouple.json")))
rows = DATA["rows"]
PASS = DATA["CMB_pass"]
TARGET = DATA["target_cluster_a2"]

plt.rcParams.update({"font.size": 13, "axes.titlesize": 13.5, "axes.labelsize": 13})


def demean(y):
    y = np.asarray(y, float)
    return y - y.mean()


def tag(ax, x, y, text, color, ha="left", va="center", fs=11.0):
    ax.text(
        x, y, text, color=color, fontsize=fs, va=va, ha=ha, zorder=6,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.88),
    )


# ---------- Figure 1: vs θ — matter ΔH/H + CMB terms (3 D panels) ----------
fig, axes = plt.subplots(3, 1, figsize=(15, 17.5), sharex=True)

for ax, r in zip(axes, rows):
    th = np.asarray(r["theta_toward_deg"])
    dHH = demean(r["dHH"])
    st = demean(r["dtt_static"])
    dop = demean(r["dtt_doppler"])
    tot = demean(r["dtt_total_static_lim"])
    Dlab = r["D_label"]
    turning = r["vrmin"] < 0

    ax.plot(th, dHH, color="C1", lw=3.4, zorder=5)
    ax.plot(th, dop, color="C0", lw=2.8, zorder=4)
    ax.plot(th, st, color="0.4", lw=2.0, ls="--", zorder=2)
    ax.plot(th, tot, color="C3", lw=1.8, ls=":", zorder=3)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 30, 60, 90, 120, 150, 180])
    ax.grid(alpha=0.25)
    ax.set_ylabel(r"anisotropy − mean")

    title = (
        "D = %s Gpc     Mc/Mcrit = %.4f     M_eff/Mcrit = %+.2e"
        % (Dlab, r["mcc"], r["Meff"])
    )
    if abs(r["static_a2"]) < 1e-5:
        title += "     endpoint static ≈ 0"
    ax.set_title(title, loc="left", pad=8)

    if Dlab == "22.28":
        ymax = 3.05
    elif Dlab == "14.7":
        ymax = 0.55
    else:
        ymax = max(np.max(np.abs(dop)), np.max(np.abs(dHH)), np.max(np.abs(st))) * 1.55
    ax.set_ylim(-ymax, ymax)

    if turning:
        ax.axvspan(0, 25, color="0.88", alpha=0.5, zorder=0)
        tag(ax, 12, 0.72 * ymax if Dlab != "22.28" else 2.35,
            "turnaround\nregion", "0.35", ha="center", fs=10)

    if Dlab == "14.7":
        tag(ax, 55, 0.38,
            "cluster ΔH/H\na₂ = %.3f (Migkas match)" % r["cluster_a2"],
            "C1", ha="left")
        tag(ax, 55, -0.22,
            "CMB Doppler\na₂ = %+.3f" % r["doppler_a2"],
            "C0", ha="left")
        tag(ax, 145, 0.12, "endpoint static\na₂ = %+.3f" % r["static_a2"],
            "0.35", ha="center")
        tag(ax, 145, -0.38, "CMB total\n(static Φ limit)", "C3", ha="center")
    elif Dlab == "22.28":
        tag(ax, 70, 2.35,
            "cluster ΔH/H\na₂ = %.3f  ← matches Migkas" % r["cluster_a2"],
            "C1", ha="left")
        tag(ax, 70, 1.15,
            "CMB Doppler still large\na₂ = %+.3f   (%.0f× over 1e−5)"
            % (r["doppler_a2"], r["ratio_doppler_a2_over_CMB_pass"]),
            "C0", ha="left")
        tag(ax, 150, 0.55, "endpoint ≈ 0\nISW ≈ 0\n(static null)",
            "0.35", ha="center")
        tag(ax, 130, -1.65, "CMB total ≈ Doppler", "C3", ha="center")
    else:
        tag(ax, 52, 0.42,
            "cluster ΔH/H\na₂ = %.3f" % r["cluster_a2"],
            "C1", ha="left")
        tag(ax, 52, -0.28,
            "CMB Doppler\na₂ = %+.3f" % r["doppler_a2"],
            "C0", ha="left")
        tag(ax, 18, 0.20, "endpoint\na₂ = %+.3f" % r["static_a2"],
            "0.35", ha="left")
        tag(ax, 150, -0.42, "CMB total", "C3", ha="center")

axes[0].plot([], [], color="C1", lw=3.4, label=r"cluster $\Delta H/H$  (matter, $R_s=0.31$ Gpc)")
axes[0].plot([], [], color="C0", lw=2.8, label=r"CMB Doppler  $-n\cdot\Delta v/c$")
axes[0].plot([], [], color="0.4", lw=2.0, ls="--", label=r"CMB endpoint  $\Phi/c^2$")
axes[0].plot([], [], color="C3", lw=1.8, ls=":", label=r"CMB total (static-$\Phi$ limit, ISW$\,=\,0$)")
axes[0].legend(loc="upper right", fontsize=10.2, framealpha=0.92)

axes[-1].set_xlabel(
    r"angle $\theta$ from CoE axis (degrees)      "
    r"$0^\circ$ = toward CoE      ·      $180^\circ$ = away from CoE"
)
fig.suptitle(
    r"Same $\Phi$: cluster $\Delta H/H$ vs CMB $\Delta T/T$  (Chris CoE compact core + $\Omega_s=0.049$ shell)"
    "\n"
    r"Nearly-static $\Phi$: matter accumulates $\Delta v=-\!\int\nabla\Phi\,dt$;  "
    r"ISW $\propto\Phi̇$ is small;  endpoint $\propto M_{\rm eff}$"
    "\n"
    r"Newtonian integrator · $\langle H\rangle=71.6$ · $R=14.3$ Gpc · Mc re-solved for cluster $a_2=0.084$",
    fontsize=14.2, y=0.995,
)
fig.tight_layout(rect=(0, 0.01, 1, 0.915))
out1 = os.path.join(HERE, "decouple_vs_theta.png")
fig.savefig(out1, dpi=130)
print("wrote", out1)


# ---------- Figure 2: |a2| bar comparison ----------
fig2, axes2 = plt.subplots(1, 3, figsize=(16.8, 6.2), sharey=False)
cats = ["ΔH/H\n(cluster)", "endpoint\nΦ/c²", "ISW\n(Φ̇≈0)", "Doppler", "CMB total\n(static Φ)"]
colors = ["C1", "0.45", "0.65", "C0", "C3"]

for ax, r in zip(axes2, rows):
    vals = [
        abs(r["cluster_a2"]),
        abs(r["static_a2"]),
        abs(r["isw_static_a2"]),
        abs(r["doppler_a2"]),
        abs(r["total_static_lim_a2"]),
    ]
    x = np.arange(len(cats))
    bars = ax.bar(x, vals, color=colors, width=0.72, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(cats, fontsize=10.5)
    ax.set_yscale("log")
    ax.set_ylim(1e-10, 2.0)
    ax.axhline(PASS, color="g", lw=2.0, ls="--", zorder=2)
    ax.axhline(TARGET, color="C1", lw=1.6, ls=":", zorder=2)
    ax.grid(alpha=0.3, which="both", axis="y")
    ax.set_title("D = %s Gpc" % r["D_label"], loc="left", pad=8)

    # in-panel labels
    tag(ax, 0.15, TARGET * 1.8, "Migkas a₂=0.084", "C1", ha="left", fs=9.5)
    tag(ax, 3.55, PASS * 2.2, "CMB limit 1e−5", "g", ha="right", fs=9.5)

    for xi, v in zip(x, vals):
        if v < 1e-9:
            lab = "≈0"
        elif v >= 0.01:
            lab = "%.3f" % v
        else:
            lab = "%.1e" % v
        ax.text(xi, max(v * 1.35, PASS * 3), lab, ha="center", va="bottom",
                fontsize=9.5, color="k", zorder=6)

axes2[0].set_ylabel(r"$|a_2|$  (log scale)")
fig2.suptitle(
    r"Decoupling scorecard: $|a_2|$ of $\Delta H/H$ vs $\Delta T/T$ pieces from the same $\Phi$"
    "\n"
    r"Static-$\Phi$ limit (ISW$\,=\,0$).  At $D=22.28$ endpoint cancels; Doppler remains $\sim 4\times 10^{4}\times$ over CMB while clusters still match.",
    fontsize=13.5,
)
fig2.tight_layout(rect=(0, 0.02, 1, 0.88))
out2 = os.path.join(HERE, "decouple_a2_bars.png")
fig2.savefig(out2, dpi=130)
print("wrote", out2)


# ---------- Figure 3: growing-mode contrast (optional, one row of 3) ----------
fig3, axes3 = plt.subplots(1, 3, figsize=(16.8, 5.8), sharey=False)
for ax, r in zip(axes3, rows):
    th = np.asarray(r["theta_toward_deg"])
    dHH = demean(r["dHH"])
    dop = demean(r["dtt_doppler"])
    isw = demean(r["dtt_isw_growing"])
    totg = demean(r["dtt_total_growing"])
    ax.plot(th, dHH, color="C1", lw=2.6)
    ax.plot(th, dop, color="C0", lw=2.4)
    ax.plot(th, isw, color="C2", lw=2.2, ls="--")
    ax.plot(th, totg, color="C3", lw=2.0, ls=":")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 60, 120, 180])
    ax.grid(alpha=0.25)
    ax.set_xlabel(r"$\theta$ (deg)")
    ax.set_title("D = %s Gpc" % r["D_label"], loc="left")
    ymax = max(0.15, float(np.max(np.abs(np.concatenate([dHH, dop, isw, totg])))) * 1.35)
    if r["D_label"] == "22.28":
        ymax = 3.05
    ax.set_ylim(-ymax, ymax)
    tag(ax, 95, 0.72 * ymax,
        "ISW grow a₂=%+.3f\nnet grow+Dop a₂=%+.3f"
        % (r["isw_growing_a2"], r["total_growing_a2"]),
        "C2", ha="center", fs=9.5)

axes3[0].set_ylabel(r"anisotropy − mean")
axes3[0].plot([], [], color="C1", lw=2.6, label=r"cluster $\Delta H/H$")
axes3[0].plot([], [], color="C0", lw=2.4, label="CMB Doppler")
axes3[0].plot([], [], color="C2", lw=2.2, ls="--", label=r"CMB ISW (growing $\Phi$, Om=0.049)")
axes3[0].plot([], [], color="C3", lw=2.0, ls=":", label="CMB Doppler + growing net")
axes3[0].legend(loc="lower right", fontsize=9.0, framealpha=0.9)
fig3.suptitle(
    r"Contrast only: growing-mode history (not Chris's static table).  "
    r"At $D=22.28$, $M_{\rm eff}\approx 0$ ⇒ SW and ISW also vanish; Doppler alone survives.",
    fontsize=12.5,
)
fig3.tight_layout(rect=(0, 0, 1, 0.90))
out3 = os.path.join(HERE, "decouple_growing_contrast.png")
fig3.savefig(out3, dpi=130)
print("wrote", out3)

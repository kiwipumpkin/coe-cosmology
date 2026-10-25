"""FULL-SIZE ΔT/T vs θ for three D. Emphasize Doppler. Labels in open space, no leader lines."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = json.load(open(os.path.join(HERE, "dtt.json")))
rows = DATA["rows"]
HRC = DATA["HRC"]

plt.rcParams.update({"font.size": 13, "axes.titlesize": 13.5, "axes.labelsize": 13})


def demean(y):
    y = np.asarray(y, float)
    return y - y.mean()


def tag(ax, x, y, text, color, ha="left", va="center", fs=11.5):
    ax.text(x, y, text, color=color, fontsize=fs, fontweight="normal",
            va=va, ha=ha, zorder=6,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.85))


fig, axes = plt.subplots(3, 1, figsize=(15, 17), sharex=True)

for ia, (ax, r) in enumerate(zip(axes, rows)):
    th = np.asarray(r["theta_toward_deg"])
    st = demean(r["dtt_static"])
    dop = demean(r["dtt_doppler"])
    tot = demean(r["dtt_total"])
    Dlab = r["D_label"]
    a2d = r["doppler_a2"]
    a2s = r["static_a2"]
    ptp = r["doppler_ptp"]
    Me = r["Meff"]
    turning = r["vrmin"] < 0

    ax.plot(th, st, color="0.4", lw=2.0, ls="--", zorder=2)
    ax.plot(th, tot, color="C3", lw=2.0, ls=":", zorder=3)
    ax.plot(th, dop, color="C0", lw=3.6, zorder=4)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 30, 60, 90, 120, 150, 180])
    ax.grid(alpha=0.25)

    title = "D = %s Gpc     Mc/Mcrit = %.4f     M_eff/Mcrit = %+.2e" % (Dlab, r["mcc"], Me)
    if abs(a2s) < 1e-5:
        title += "     static ≈ 0 (exact cancel)"
    ax.set_title(title, loc="left", pad=8)

    if Dlab == "22.28":
        ymax = 3.05
    elif Dlab == "14.7":
        ymax = 0.55
    else:
        ymax = max(np.max(np.abs(dop)), np.max(np.abs(st)), np.max(np.abs(tot))) * 1.45
    ax.set_ylim(-ymax, ymax)
    ax.set_ylabel(r"$\Delta T/T$ − mean")

    if turning:
        ax.axvspan(0, 25, color="0.88", alpha=0.5, zorder=0)
        # sit ABOVE zero in the shaded band, clear of the spike
        tag(ax, 12, 0.78 * ymax if Dlab != "22.28" else 2.4,
            "turnaround\nregion", "0.35", ha="center", fs=10)

    if Dlab == "14.7":
        tag(ax, 58, 0.40,
            "Doppler (emitter motion)\npeak-to-peak = %.3f\na₂ = %+.3f" % (ptp, a2d),
            "C0", ha="left")
        tag(ax, 148, 0.10, "static\na₂ = %+.3f" % a2s, "0.35", ha="center")
        tag(ax, 148, -0.32, "total", "C3", ha="center")
    elif Dlab == "22.28":
        tag(ax, 75, 1.35,
            "Doppler (emitter motion)\npeak-to-peak = %.3f\na₂ = %+.3f" % (ptp, a2d),
            "C0", ha="left")
        tag(ax, 155, 0.55, "static ≈ 0\n(core = cavity)", "0.35", ha="center")
        tag(ax, 130, -1.55, "total ≈ Doppler", "C3", ha="center")
    else:
        tag(ax, 58, 0.36,
            "Doppler (emitter motion)\npeak-to-peak = %.3f\na₂ = %+.3f" % (ptp, a2d),
            "C0", ha="left")
        tag(ax, 18, 0.18, "static\na₂ = %+.3f" % a2s, "0.35", ha="left")
        tag(ax, 155, -0.40, "total", "C3", ha="center")

axes[0].plot([], [], color="C0", lw=3.6, label="Doppler from emitter motion  (target residual)")
axes[0].plot([], [], color="0.4", lw=2.0, ls="--", label="static gravitational redshift")
axes[0].plot([], [], color="C3", lw=2.0, ls=":", label="total = static + Doppler")
axes[0].legend(loc="upper right", fontsize=10.5, framealpha=0.92)

axes[-1].set_xlabel(
    r"angle $\theta$ from CoE axis (degrees)      "
    r"$0^\circ$ = toward CoE      ·      $180^\circ$ = away from CoE"
)

fig.suptitle(
    "CMB emission sphere (R = 14.3 Gpc):  ΔT/T vs angle\n"
    "Chris CoE compact core + Ωₛ = 0.049 shell  ·  Newtonian, today's field  ·  ⟨H⟩ = 71.6\n"
    "Thick blue = kinematic Doppler residual he wants to eliminate with new physics / interpretation",
    fontsize=14.5, y=0.995,
)
fig.tight_layout(rect=(0, 0.01, 1, 0.925))
out = os.path.join(HERE, "dtt_vs_theta.png")
fig.savefig(out, dpi=130)
print("wrote", out)

# Fractional caveat: 3 panels
fig2, axes2 = plt.subplots(1, 3, figsize=(16.5, 5.8), sharey=False)
for ax, r in zip(axes2, rows):
    th = np.asarray(r["theta_toward_deg"])
    dop = demean(r["dtt_doppler"])
    frac = demean(r["dtt_doppler_frac"])
    ax.plot(th, dop, color="C0", lw=2.8)
    ax.plot(th, frac, color="C2", lw=2.4, ls="--")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 60, 120, 180])
    ax.grid(alpha=0.25)
    ax.set_xlabel(r"$\theta$ (deg)")
    ymax = 3.05 if r["D_label"] == "22.28" else max(np.max(np.abs(dop)), np.max(np.abs(frac))) * 1.4
    ax.set_ylim(-ymax, ymax)
    ax.set_title("D = %s Gpc" % r["D_label"], loc="left")
    # labels in open space
    if r["D_label"] == "22.28":
        tag(ax, 95, 2.2, "linear\na₂=%+.3f\nptp=%.2f" % (r["doppler_a2"], r["doppler_ptp"]), "C0", ha="center", fs=10.5)
        tag(ax, 95, -1.8, "fractional\na₂=%+.3f\nptp=%.2f" % (r["frac_a2"], r["frac_ptp"]), "C2", ha="center", fs=10.5)
    else:
        tag(ax, 100, 0.70 * ymax, "linear\na₂=%+.3f" % r["doppler_a2"], "C0", ha="center", fs=10.5)
        tag(ax, 100, -0.55 * ymax, "fractional\na₂=%+.3f" % r["frac_a2"], "C2", ha="center", fs=10.5)
axes2[0].set_ylabel(r"$\Delta T/T$ − mean")
axes2[0].plot([], [], color="C0", lw=2.8, label="linear Doppler −n·Δv/c")
axes2[0].plot([], [], color="C2", lw=2.4, ls="--", label="fractional = linear / (HR/c)")
axes2[0].legend(loc="lower right", fontsize=9.5, framealpha=0.9)
fig2.suptitle(
    "Caveat: HR/c = %.2f  ⇒  linear v/c outside its range.  "
    "Fractional anisotropic-expansion estimate still has |a₂| ≫ 1e−5 (does not save the model)." % HRC,
    fontsize=12.5,
)
fig2.tight_layout(rect=(0, 0, 1, 0.90))
out2 = os.path.join(HERE, "dtt_doppler_linear_vs_frac.png")
fig2.savefig(out2, dpi=130)
print("wrote", out2)

"""Residual panel + percentages for the 4-clump tidal fit.

Reuses saved clump-fit.json parameters and cmb-cancel-story.json run-2 params.
Does NOT refit.
"""
import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

src = open("cmb-cancel.py").read().split("if __name__")[0]
g = {}; exec(src, g)
story = json.load(open("cmb-cancel-story.json"))
P = story["params"]; phi = P[4]
ELL = g["ELL"]; G = g["G_FULL"]

R0 = G - 1.0  # starting shell-geometry anisotropy (before SR + mass)
R = g["resid"](G, g["S_model"](*P)[0]); R = R - R.mean()  # after SR + mass (run 2)

cf = json.load(open("clump-fit.json"))
p = np.array(cf["2 pairs"]["p"])  # reuse saved 2-pair / 4-clump params

th = np.radians(ELL)


def field(alpha, D):
    x, y = np.cos(th), np.sin(th)
    X, Y = D * np.cos(alpha), D * np.sin(alpha)
    dx, dy = X - x, Y - y
    r3 = (dx * dx + dy * dy) ** 1.5
    ar = (dx * x + dy * y) / r3
    a0x, a0y = X / D ** 3, Y / D ** 3
    ar = ar - (a0x * x + a0y * y)
    return ar - ar.mean()


def model_2pairs(p):
    D = p[0]; out = np.zeros_like(th); i = 1
    for _ in range(2):
        k_, off = p[i], p[i + 1]; i += 2
        out += k_ * (field(np.radians(phi + off), D) + field(np.radians(phi - off), D))
    return out


m = model_2pairs(p)
left = R - m


def rms(x):
    return float(np.sqrt(np.mean(x ** 2)))


a, b, c = rms(R0), rms(R), rms(left)
max_abs_left = float(np.abs(left).max())

sr_amp = 100 * (a - b) / a
sr_var = 100 * (a ** 2 - b ** 2) / a ** 2
tot_amp = 100 * (a - c) / a
tot_var = 100 * (a ** 2 - c ** 2) / a ** 2
clump_of_start_amp = 100 * (b - c) / a
clump_of_start_var = 100 * (b ** 2 - c ** 2) / a ** 2
clump_of_post_amp = 100 * (b - c) / b
clump_of_post_var = 100 * (1 - c ** 2 / b ** 2)

rad = np.radians(ELL)
var_left = float(np.mean(left ** 2))
modes = []
for mm in range(1, 13):
    am = float(2 * np.mean(left * np.cos(mm * rad)))
    bm = float(2 * np.mean(left * np.sin(mm * rad)))
    Am = float(np.hypot(am, bm))
    vf = float(0.5 * Am ** 2 / var_left) if var_left > 0 else 0.0
    modes.append(dict(m=mm, A=Am, var_frac=vf))

D = float(p[0])
angs = [(phi + p[2]) % 360, (phi - p[2]) % 360, (phi + p[4]) % 360, (phi - p[4]) % 360]
angs_sorted = sorted(angs)

out = dict(
    rms_starting_geometry=a,
    rms_after_SR_mass=b,
    rms_after_4_clumps=c,
    max_abs_residual_after_4_clumps=max_abs_left,
    percent_of_starting=dict(
        SR_mass_amplitude=sr_amp,
        SR_mass_variance=sr_var,
        clumps_alone_amplitude=clump_of_start_amp,
        clumps_alone_variance=clump_of_start_var,
        total_after_clumps_amplitude=tot_amp,
        total_after_clumps_variance=tot_var,
    ),
    percent_of_post_SR_residual=dict(
        clumps_amplitude=clump_of_post_amp,
        clumps_variance=clump_of_post_var,
    ),
    clump_params_reused=dict(
        p=p.tolist(), D_in_visible_radii=D, D_Gpc=D * 14.3,
        angles_deg=angs_sorted, mass_kg_each=7.7e51,
    ),
    leftover_fourier_modes=modes,
)
json.dump(out, open("clump-fit-residual.json", "w"), indent=1)

plt.rcParams.update({"font.size": 15})
fig = plt.figure(figsize=(18, 15))
gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 1.05],
              hspace=0.72, wspace=0.30,
              top=0.92, bottom=0.06, left=0.07, right=0.98)

ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(ELL, R * 1e3, "C3", lw=2.5, label="residual after special relativity (SR) + mass")
ax1.plot(ELL, m * 1e3, "k--", lw=2, label="4 clumps, tidal pull")
ax1.set_xlim(0, 360); ax1.grid(alpha=0.3)
ax1.set_xlabel(r"Galactic longitude $\ell$ (deg)")
ax1.set_ylabel(r"$\Delta T/T\times10^{3}$")
ax1.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1, frameon=False, fontsize=14)
ax1.set_title(
    f"SR + mass removes {sr_amp:.1f}% of starting residual (RMS); "
    f"4 clumps remove a further {clump_of_start_amp:.1f}%; total {tot_amp:.1f}%\n"
    f"by variance: {sr_var:.2f}% → total {tot_var:.3f}%; "
    f"4 clumps explain {clump_of_post_var:.0f}% of the post-SR residual variance",
    pad=72, fontsize=14,
)

ax2 = fig.add_subplot(gs[0, 1], projection="polar")
ax2.plot(np.linspace(0, 2 * np.pi, 400), np.ones(400), "C0", lw=2,
         label="visible universe, 14.3 Gpc")
ax2.scatter(np.radians(angs), [D] * 4, s=400, color="C3", zorder=3,
            label=f"clumps at {D * 14.3:.1f} Gpc, 7.7e51 kg each")
ax2.plot([np.radians(phi)] * 2, [0, 2.1], "k", lw=1.5)
ax2.plot([np.radians(phi + 180)] * 2, [0, 2.1], "k--", lw=1.5)
for A in angs:
    ax2.text(np.radians(A), D + 0.32, f"{A:.0f}°", ha="center", va="center", fontsize=14)
ax2.set_ylim(0, 2.35); ax2.set_yticklabels([])
ax2.legend(loc="lower center", bbox_to_anchor=(0.5, 1.06), frameon=False, fontsize=14)
ax2.set_title("4-clump placement about the 15.16° axis", pad=95, fontsize=15)

ax3 = fig.add_subplot(gs[1, :])
ax3.plot(ELL, left * 1e3, "C0", lw=2.4, label="leftover after 4 clumps (residual − model)")
ax3.axhline(0, color="k", lw=0.7)
cmb_scaled = 1e-5 * 1e3  # 0.01 on this axis — sits on zero vs leftover ~±2
ax3.axhline(cmb_scaled, ls="--", color="purple", lw=1.8)
ax3.axhline(-cmb_scaled, ls="--", color="purple", lw=1.8,
            label=r"CMB level $\pm10^{-5}$ (near zero on this scale)")
ax3.set_xlim(0, 360); ax3.grid(alpha=0.3)
ax3.set_xlabel(r"Galactic longitude $\ell$ (deg)")
ax3.set_ylabel(r"$\Delta T/T\times10^{3}$")
ymin, ymax = left.min() * 1e3, left.max() * 1e3
pad = 0.18 * max(abs(ymin), abs(ymax))
ax3.set_ylim(ymin - pad, ymax + pad)
ax3.legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=2, frameon=False, fontsize=14)
ax3.set_title(
    f"Residual after 4 clumps: RMS = {c:.2e}  (max |ΔT/T| = {max_abs_left:.2e})\n"
    f"Cumulative from start: {tot_amp:.1f}% of starting RMS removed "
    f"({tot_var:.3f}% of starting variance); "
    f"{clump_of_post_var:.0f}% of post-SR variance removed by the clumps",
    pad=52, fontsize=14,
)

fig.savefig("clump-fit-residual.png", dpi=110, bbox_inches="tight")
print("a,b,c =", a, b, c)
print("max|left| =", max_abs_left)
print("saved clump-fit-residual.png / .json")

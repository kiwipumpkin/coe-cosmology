#!/usr/bin/env python3
"""Standalone Amp x phi R^2 heat maps (A and B) for the ST+CoE paper.

A: dense grid recomputed. Baseline A_phys recovered exactly from
   best_fit_2ell_multiplier.npz: A_phys = A_before / (1 - Amp1*sin(l+phi1)).
   R^2_A = affine least-squares R^2 of (1 - Amp*sin(l+phi))*A_phys vs Migkas A/A_all
   (71 digitised points). Validated against all 35 stored level-set grid points in
   rotation_Amp_phi_R2_ranked.csv (max |dR^2| printed below).
B: plotted directly from the stored 13 x 37 grid rotation_Amp_phi_B_R2_ranked.csv
   (no recomputation).
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = "/workspace/cosmology/map_compare/amp_mean_param_noise"
OUT = "/workspace/cosmology/figures"
plt.rcParams.update({"font.size": 11})

# ---------------- A ----------------
d = np.load(f"{SRC}/best_fit_2ell_multiplier.npz")
lm, A = d["ell_m"], d["A_m"]
Aphys = d["A_before"] / (1 - d["Amp1"] * np.sin(np.radians(lm) + d["phi1"]))
sst = ((A - A.mean()) ** 2).sum()

def r2A(amp, phi_rad):
    m = (1 - amp * np.sin(np.radians(lm) + phi_rad)) * Aphys
    X = np.column_stack([m, np.ones_like(m)])
    c, *_ = np.linalg.lstsq(X, A, rcond=None)
    return 1 - ((A - X @ c) ** 2).sum() / sst

stored = np.genfromtxt(f"{SRC}/rotation_Amp_phi_R2_ranked.csv", delimiter=",", names=True)
maxdiff = max(abs(r["R2"] - r2A(r["Amp"], r["phi_rad"])) for r in stored)

amps = np.round(np.arange(0.0, -0.1401, -0.01), 3)
phis = np.arange(0.0, 180.01, 5.0)
Z = np.array([[r2A(a, np.radians(p)) for p in phis] for a in amps])
ib = np.unravel_index(np.argmax(Z), Z.shape)

def edges(v):
    dv = np.diff(v).mean()
    return np.concatenate([v - dv / 2, [v[-1] + dv / 2]])

fig, ax = plt.subplots(figsize=(6.4, 4.4))
pm = ax.pcolormesh(edges(phis), edges(amps), Z, cmap="magma", shading="flat", edgecolors="face", linewidth=0, rasterized=True)
ax.set_ylim(0.005, -0.145)
ax.plot(20.05, -0.06, marker="^", ms=10, color="#17becf", mec="k", mew=0.6, ls="none",
        label=r"Amp $=-0.06$, $\varphi=20^\circ$")
ax.plot(phis[ib[1]], amps[ib[0]], marker="x", ms=9, mew=2, color="k", ls="none",
        label=rf"grid max: Amp $={amps[ib[0]]:.2f}$, $\varphi={phis[ib[1]]:.0f}^\circ$, $R^2={Z[ib]:.3f}$")
ax.set_xlabel(r"rotation phase $\varphi$ [deg]")
ax.set_ylabel(r"Amp")
cb = fig.colorbar(pm, ax=ax)
cb.set_label(r"$R^2_A$ vs Migkas $A/A_{\mathrm{all}}$")
ax.legend(loc="upper right", fontsize=8.5, framealpha=0.9)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/heatmap_R2_A_Amp_phi.{ext}", dpi=300)
plt.close(fig)

# ---------------- B ----------------
b = np.genfromtxt(f"{SRC}/rotation_Amp_phi_B_R2_ranked.csv", delimiter=",", names=True)
bamps = np.array(sorted(set(np.round(b["Amp_B"], 4))))
bphis = np.array(sorted(set(np.round(b["phi_deg"], 2))))
ZB = np.full((len(bamps), len(bphis)), np.nan)
for r in b:
    ZB[np.searchsorted(bamps, round(r["Amp_B"], 4)), np.searchsorted(bphis, round(r["phi_deg"], 2))] = r["R2"]
assert not np.isnan(ZB).any()
best = b[np.argmax(b["R2"])]

fig, ax = plt.subplots(figsize=(6.4, 4.4))
pm = ax.pcolormesh(edges(bphis), edges(bamps), ZB, cmap="plasma", shading="flat", edgecolors="face", linewidth=0, rasterized=True)
ax.axhline(0.06, color="#7CFC00", ls=":", lw=1.4, label=r"Amp$_B=+0.06$ ($=-$Amp$_A$)")
ax.axvline(11.46, color="#17becf", ls=":", lw=1.2, label=r"$\varphi=11.5^\circ$ (best for $A$)")
ax.plot(best["phi_deg"], best["Amp_B"], marker="*", ms=15, color="w", mec="k", mew=0.8, ls="none",
        label=rf"best: Amp$_B={best['Amp_B']:.2f}$, $\varphi={best['phi_deg']:.0f}^\circ$, $R^2={best['R2']:.3f}$")
ax.set_xlabel(r"rotation phase $\varphi$ [deg]")
ax.set_ylabel(r"Amp$_B$")
cb = fig.colorbar(pm, ax=ax)
cb.set_label(r"$R^2_B$ vs Migkas $B/B_{\mathrm{all}}$")
ax.legend(loc="lower right", fontsize=8.5, framealpha=0.9)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/heatmap_R2_B_AmpB_phi.{ext}", dpi=300)
plt.close(fig)

meta = {
    "A": {"method": "recomputed dense grid from npz baseline; affine R^2",
          "amp_grid": amps.tolist(), "phi_grid_deg": phis.tolist(),
          "validation_max_abs_dR2_vs_stored_35pts": maxdiff,
          "grid_max": {"Amp": float(amps[ib[0]]), "phi_deg": float(phis[ib[1]]), "R2": float(Z[ib])},
          "R2_at_Amp-0.06_phi0.29rad": r2A(-0.06, 0.29),
          "R2_at_Amp-0.06_phi11.46deg": r2A(-0.06, np.radians(11.46))},
    "B": {"method": "plotted directly from stored CSV grid", "n": int(len(b)),
          "best": {"Amp_B": float(best["Amp_B"]), "phi_deg": float(best["phi_deg"]), "R2": float(best["R2"])}},
}
json.dump(meta, open(f"{OUT}/heatmaps_R2_Amp_phi.json", "w"), indent=2)
print(json.dumps({k: {kk: vv for kk, vv in v.items() if "grid" not in kk or kk == "grid_max"} for k, v in meta.items()}, indent=1))

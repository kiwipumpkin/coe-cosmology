"""Pure exponential mass density vs distance from CoE (Centre of Expansion).

rho(r) = rho_0 * exp(-r/L), with M = M(<D) setting the scale:
  f(u) = 1 - exp(-u)*(1 + u + u^2/2)
  rho_0 = M / (8 pi L^3 f(D/L))

Plotted: rho(r)/rho(D) = exp(-(r-D)/L)  (M cancels; shape set by L only).
No ptshell shell edge / cliff.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

OUT_PNG = Path("/workspace/cosmology-course/notes/exponential-density-vs-D.png")
OUT_CSV = Path("/workspace/cosmology-course/notes/exponential-density-vs-D.csv")

# Reference Earth–CoE distance (ptshell WLS D, for apples-to-apples marking)
D = 58.2  # Gpc
R_NEAR = 3.0  # Gpc (near-zone / CatWISE-ish scale)
R_VIS = 14.3  # Gpc (vision / Hubble-sphere scale used in Migkas geometry)

# L/D = 0.17 ≈ 1/5.9 (no-boxcar / extrema baseline)
# Plus L that give pure-exp contrast delta = 2 sinh(R/L) ≈ 0.01–0.03
def L_for_delta(R, delta):
    return R / np.arcsinh(delta / 2.0)


LD_BASE = 0.17
L_curves = [
    (LD_BASE * D, f"L/D = {LD_BASE:.2f} (no-boxcar baseline)", "C0", "-"),
    (0.18 * D, "L/D = 0.18 (spreadsheet D/R=5 geometry)", "C1", "--"),
    (L_for_delta(R_NEAR, 0.02), rf"$\delta(R\!=\!{R_NEAR})=2\sinh(R/L)=0.02$", "C2", "-"),
    (L_for_delta(R_VIS, 0.02), rf"$\delta(R\!=\!{R_VIS})=2\sinh(R/L)=0.02$", "C3", "-."),
    (L_for_delta(R_NEAR, 0.01), rf"$\delta(R\!=\!{R_NEAR})=0.01$", "C4", ":"),
    (L_for_delta(R_NEAR, 0.03), rf"$\delta(R\!=\!{R_NEAR})=0.03$", "C5", (0, (3, 1, 1, 1))),
]

r = np.linspace(0.0, D + 1.6 * R_VIS, 2500)  # from CoE outward
rows = ["r_Gpc," + ",".join(f"rho_over_rhoD_L={L:.6g}" for L, *_ in L_curves)]

fig, ax = plt.subplots(figsize=(10.5, 6.2))
for L, lab, color, ls in L_curves:
    y = np.exp(-(r - D) / L)  # rho(r)/rho(D)
    ax.plot(r, y, color=color, ls=ls, lw=2.2, label=rf"{lab}; $L={L:.1f}$ Gpc")
    rows.append(",".join([f"{r[0]:.6g}"] + []) )  # placeholder; fill below

# CSV: full columns
csv_data = [r] + [np.exp(-(r - D) / L) for L, *_ in L_curves]
header = "r_Gpc," + ",".join(f"rho_over_rhoD_L{L:.4g}" for L, *_ in L_curves)
np.savetxt(OUT_CSV, np.column_stack(csv_data), delimiter=",", header=header, comments="")

ax.axvline(D, color="k", lw=1.6, ls="--", label=rf"Earth at $D={D:.1f}$ Gpc (CoE→Earth)")
ax.axvspan(D - R_VIS, D + R_VIS, color="C2", alpha=0.12, label=rf"Vision sphere $D\pm R$, $R={R_VIS}$ Gpc")
ax.axvline(D - R_NEAR, color="0.35", lw=1.0, ls=":")
ax.axvline(D + R_NEAR, color="0.35", lw=1.0, ls=":", label=rf"$D\pm {R_NEAR:g}$ Gpc")
ax.axhline(1.0, color="0.5", lw=0.8, alpha=0.7)

ax.set_xlabel(r"Distance from CoE (Centre of Expansion), $r$ (Gpc)", fontsize=13)
ax.set_ylabel(r"Mass density ratio $\rho(r)/\rho(D)$", fontsize=13)
ax.set_title(
    r"Pure exponential CoE mass: $\rho\propto e^{-r/L}$ (no shell cliff)"
    "\n"
    rf"Shape set by $L$ only; $M\equiv M(<D)$ cancels in $\rho/\rho(D)$.  "
    rf"Contrast $\delta(R)=2\sinh(R/L)$.",
    fontsize=14,
)
ax.legend(fontsize=14, loc="upper right", framealpha=0.92)
ax.set_xlim(0, D + 1.55 * R_VIS)
ax.set_ylim(0, 3.2)
ax.tick_params(labelsize=12)
ax.grid(alpha=0.25)
# Annotate CatWISE-scale note
ax.text(
    0.02,
    0.02,
    "CatWISE-scale $\\delta\\sim1$–3% needs huge $L$ (hundreds–thousands of Gpc)\n"
    "or a uniform floor that dilutes $f_{\\mathrm{CoE}}$; baseline $L/D\\approx0.17$ gives large $\\delta$.",
    transform=ax.transAxes,
    fontsize=11,
    va="bottom",
    ha="left",
    bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.85),
)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=160)
print("wrote", OUT_PNG)
print("wrote", OUT_CSV)
for L, lab, *_ in L_curves:
    d3 = 2 * np.sinh(R_NEAR / L)
    d14 = 2 * np.sinh(R_VIS / L)
    print(f"L={L:.3f} Gpc  LD={L/D:.4f}  delta(3)={d3:.4g}  delta(14.3)={d14:.4g}  | {lab}")

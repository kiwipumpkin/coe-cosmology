#!/usr/bin/env python3
import json, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "/workspace/cosmology/figures"
J = json.load(open(f"{OUT}/R2_vs_Dr_1to10.json")); R = J["results"]
fig, ax = plt.subplots(figsize=(8, 5))
for fam, strict, lab, c, ls, mk in [("a_exp_floor", False, r"exp core + floor ($\kappa$ free)", "#0072B2", "-", "o"),
                                    ("a_exp_floor", True, r"exp core + floor, $\rho(D)=0.3\rho_c$ fixed", "#D55E00", "--", "s"),
                                    ("d_pure_exp", False, r"pure exponential ($\kappa$ free)", "#CC79A7", ":", "^")]:
    rr = sorted([r for r in R if r["family"] == fam and r["strict_local_density"] == strict], key=lambda r: r["D_over_r"])
    ax.plot([r["D_over_r"] for r in rr], [r["R2_direct"] for r in rr], ls=ls, marker=mk, color=c, ms=5, label=lab)
ax.axvspan(3, 4, color="0.88", zorder=0, label="D/r = 3–4"); ax.axhline(0.75, color="0.5", ls="--", lw=0.8)
ax.set_xscale("log"); ax.set_xticks([1.05, 1.2, 1.5, 2, 3, 4, 5, 7, 10]); ax.set_xticklabels(["1.05", "1.2", "1.5", "2", "3", "4", "5", "7", "10"])
ax.set_xlabel(r"$D/r$  ($d=r$, $d/D = r/D$)"); ax.set_ylabel(r"best $R^2_A$ (direct, no rescale) vs Migkas Fig. 6")
top = ax.secondary_xaxis("top", functions=(lambda x: 14.3 * x, lambda x: x / 14.3)); top.set_xlabel(r"$D$ [Gpc]  ($r=14.3$ Gpc)"); top.set_xticks([15, 20, 30, 40, 60, 80, 100, 143]); top.set_xticklabels(["15", "20", "30", "40", "60", "80", "100", "143"]); top.minorticks_off()
ax.set_ylim(0.6, 0.8)
rs = [r for r in R if r["family"] == "a_exp_floor" and r["strict_local_density"] and r["D_over_r"] == 1.05][0]
ax.annotate(f"strict, D/r=1.05: R²={rs['R2_direct']:.3f} (off scale)", xy=(1.07, 0.602), xytext=(1.12, 0.615), fontsize=7.5, color="#D55E00", arrowprops=dict(arrowstyle="->", color="#D55E00"))
ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8, loc="lower right")
ax.set_title(r"Full re-optimisation of $\kappa,\ell_0,b_0$ (|b$_0$|≤30°), $D/L$, $f_u$", fontsize=10)
fig.text(0.01, 0.005, "71 overlapping 90° sectors (N_eff≈16); R² differences are indicative.", fontsize=7, color="#444")
fig.tight_layout(rect=(0, 0.02, 1, 1))
for e in ("png", "pdf"): fig.savefig(f"{OUT}/R2_vs_Dr_1to10.{e}", dpi=200)
for r in sorted(R, key=lambda r: (r["family"], r["strict_local_density"], r["D_over_r"])):
    h = r["harmonics"]
    print(f"{r['family']:12s} {int(r['strict_local_density'])} {r['D_over_r']:5.2f} {r['D_Gpc']:6.1f} {r['D_Gly']:6.0f} R2={r['R2_direct']:.4f} "
          f"L={r['L_Gpc']:.3g} D/L={r['density']['D_over_L']:.3g} k={r['kappa']:.3f} fu={r['f_u']:.3f} b0={r['b0_deg']:.1f} l0={r['l0_deg']:.1f} "
          f"rhoD={r['rho_D_over_rhocrit']:.3g} a1={h['lsq71_a1']:.3f} a2={h['lsq71_a2']:.3f} a3={h['lsq71_a3']:.3f}")

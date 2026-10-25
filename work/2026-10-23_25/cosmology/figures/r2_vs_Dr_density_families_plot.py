#!/usr/bin/env python3
"""Plot + summary JSON/table for r2_vs_Dr_density_families.py results."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/workspace/cosmology/figures"
raw = json.load(open(f"{OUT}/R2_vs_Dr_density_families_raw.json"))
res = raw["results"]
LAB = {"a_exp_floor": "(a) exp core + floor", "b_powerlaw_floor": r"(b) $(1+R/R_c)^{-n}$ + floor",
       "c_gauss_floor": "(c) Gauss core + floor", "d_pure_exp": "(d) pure exponential"}
COL = {"a_exp_floor": "#0072B2", "b_powerlaw_floor": "#D55E00", "c_gauss_floor": "#009E73", "d_pure_exp": "#CC79A7"}

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.0), gridspec_kw=dict(width_ratios=[1.35, 1]))
ax = axes[0]
for strict, ls, mk in ((False, "-", "o"), (True, ":", "x")):
    for fam in LAB:
        rr = sorted([r for r in res if r["family"] == fam and r["strict_local_density"] == strict], key=lambda r: r["D_over_r"])
        if not rr: continue
        ax.plot([r["D_over_r"] for r in rr], [r["R2_direct"] for r in rr], ls=ls, marker=mk, color=COL[fam], ms=5,
                lw=2 if not strict else 1.3, label=LAB[fam] + ("" if not strict else r", $\kappa$ fixed by $\rho(D)=0.3\rho_c$"))
ax.axvspan(3, 4, color="0.85", zorder=0, label="D/r = 3–4")
ax.axhline(0.75, color="0.4", ls="--", lw=1)
ax.axhline(0.776, color="0.4", ls="-.", lw=0.8)
ax.text(5.02, 0.776, "0.776 ref", fontsize=7, va="center")
ax.set_xlabel(r"$D/r$  ($d = r$, so $d/D = r/D$)")
ax.set_ylabel(r"best $R^2_A$ (direct, no rescale) vs Migkas Fig. 6")
ax.set_xlim(1.1, 5.4)
top = ax.secondary_xaxis("top", functions=(lambda x: x * 14.3, lambda y: y / 14.3))
top.set_xlabel(r"$D$ [Gpc]  ($r = 14.3$ Gpc)")
ax.grid(alpha=0.3)
ax.legend(fontsize=7.2, loc="lower left", framealpha=0.92)
ax.set_title(r"Full re-optimisation of $\kappa$, $\ell_0$, $b_0$ (|b$_0$|≤30°), density params", fontsize=10)

ax = axes[1]
for fam in LAB:
    rr = sorted([r for r in res if r["family"] == fam and not r["strict_local_density"]], key=lambda r: r["D_over_r"])
    ax.plot([r["D_over_r"] for r in rr], [r["kappa"] for r in rr], "-o", ms=4, color=COL[fam], label=LAB[fam] + r": $\kappa$")
    ax.plot([r["D_over_r"] for r in rr], [r["rho_D_over_rhocrit"] for r in rr], "--s", ms=3, color=COL[fam], alpha=0.7)
ax.plot([], [], "k--s", ms=3, label=r"implied $\rho(D)/\rho_c = 0.3\,\kappa/\kappa_{\rm base}$")
ax.axhline(0.15, color="0.5", lw=0.8, ls=":"); ax.axhline(0.3, color="0.5", lw=0.8, ls=":")
ax.axvspan(3, 4, color="0.85", zorder=0)
ax.set_yscale("log"); ax.set_ylim(1e-3, 2); ax.set_xlabel(r"$D/r$"); ax.set_ylabel("best-fit value (clipped at 1e-3)")
ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=7, loc="best")
ax.set_title(r"Best-fit $\kappa$ and implied local density", fontsize=10)
fig.text(0.01, 0.005, "Neighbouring 90° Migkas sectors overlap (71 points not independent; N_eff≈16); R² differences are indicative.",
         fontsize=7, color="#444")
fig.tight_layout(rect=(0, 0.02, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/R2_vs_Dr_density_families.{ext}", dpi=200)

# table
lines = []
hdr = f"{'family':18s} {'D/r':>5s} {'d/D':>5s} {'R2':>6s} {'corr2':>6s} {'kappa':>7s} {'l0':>6s} {'b0':>6s}  density-params  rhoD/rc rho0/rhoD  a1  l1  a2  l2  a3  l3"
lines.append(hdr)
for strict in (False, True):
    lines.append(f"--- strict local density (kappa = kappa_base) = {strict} ---")
    for fam in LAB:
        for r in sorted([r for r in res if r["family"] == fam and r["strict_local_density"] == strict], key=lambda r: r["D_over_r"]):
            dp = " ".join(f"{k}={v:.3g}" for k, v in r["density"].items())
            h = r["harmonics"]
            lines.append(f"{fam:18s} {r['D_over_r']:5.2f} {r['d_over_D']:5.3f} {r['R2_direct']:6.3f} {r['R2_corr2']:6.3f} {r['kappa']:7.3f} "
                         f"{r['l0_deg']:6.1f} {r['b0_deg']:6.1f}  {dp:34s} {r['rho_D_over_rhocrit']:7.3f} {r['rho_center_over_rho_D']:9.3g} "
                         f"{h['lsq71_a1']:.3f} {h['lsq71_l1_deg']:5.1f} {h['lsq71_a2']:.3f} {h['lsq71_l2_deg']:5.1f} {h['lsq71_a3']:.3f} {h['lsq71_l3_deg']:5.1f}  nopt={r['n_local_optima_within_0.01']}")
open(f"{OUT}/R2_vs_Dr_density_families_table.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
summary = dict(doc=raw["doc"], sanity=raw["sanity"],
               data_harmonics_lsq71={"a1": 0.0530, "l1_deg": 54.9, "a2": 0.0741, "l2_deg": 14.8, "a3": 0.0320, "l3_deg": 67.0,
                                     "note": "from dq_vs_model_match.json (a0+m1+m2 fit for a1,a2; a0+m1..3 fit for a3)"},
               results=res)
json.dump(summary, open(f"{OUT}/R2_vs_Dr_density_families.json", "w"), indent=1)

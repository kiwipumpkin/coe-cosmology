"""Percent of C_l power per l, l=2..30: blended A/A_all sky (monopole removed) vs Planck 2018 TT.
Shares are of C_l (C_l = 2πD_l/(l(l+1)) for Planck), same convention as earlier share bars.
Model l=1 share of l=1..30 annotated. Reads/extends A_sky_Cl_overlay_pct.json (run A_sky_Cl_overlay_planckTT.py first)."""
import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
J = json.load(open("A_sky_Cl_overlay_pct.json"))
cl = np.array(J["Cl_model"]); lp = np.array(J["planck_l"]); Dp = np.array(J["planck_Dl_uK2"])
Cp = 2*np.pi*Dp/(lp*(lp+1)); L = np.arange(2, 31)
pm = 100*cl[2:31]/cl[2:31].sum(); pp = 100*Cp/Cp.sum()
f = 2*np.pi/(lp*(lp+1)); elo = 100*np.array(J["planck_err_lo"])*f/Cp.sum(); ehi = 100*np.array(J["planck_err_hi"])*f/Cp.sum()
s1 = 100*cl[1]/cl[1:31].sum()
fig, ax = plt.subplots(figsize=(13, 5)); w = 0.4
ax.bar(L - w/2, pm, w, color="C3", label="model (blended A/A_all sky)")
ax.bar(L + w/2, pp, w, color="k", alpha=.75, yerr=[elo, ehi], capsize=1.5, label="Planck 2018 TT")
ax.set(xlabel="ℓ", ylabel="% of C_ℓ power over ℓ=2–30", title="Percent of power per ℓ (monopole removed; dipole shown separately)")
ax.text(0.55, 0.85, f"model ℓ=1 (dipole) = {s1:.1f}% of model C_ℓ power over ℓ=1–30\n(Planck TT has no ℓ=1 entry)", transform=ax.transAxes, fontsize=9,
        bbox=dict(fc="w", ec="0.6"))
ax.legend(loc="center right", fontsize=8); ax.grid(axis="y", alpha=.3)
ins = ax.inset_axes([0.33, 0.35, 0.2, 0.5]); m = L <= 10
ins.bar(L[m] - w/2, pm[m], w, color="C3"); ins.bar(L[m] + w/2, pp[m], w, color="k", alpha=.75)
ins.set_yscale("log"); ins.set_title("ℓ=2–10, log y", fontsize=8); ins.tick_params(labelsize=7)
fig.tight_layout(); fig.savefig("A_sky_Cl_pctpower_vs_planck.png", dpi=150)
J.update({"pct_l": L.tolist(), "pct_model_l2_30": pm.tolist(), "pct_planck_l2_30": pp.tolist(), "model_l1_share_of_l1_30_pct": s1})
json.dump(J, open("A_sky_Cl_overlay_pct.json", "w"), indent=1)
print("l1 share", s1)
for i in range(9): print(L[i], round(pm[i], 2), round(pp[i], 2))

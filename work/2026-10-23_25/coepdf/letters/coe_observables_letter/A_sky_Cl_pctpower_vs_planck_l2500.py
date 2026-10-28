"""Percent of C_l power per l over l=2..2500 (log y) and cumulative %, blended A/A_all sky (nside 1024) vs Planck 2018 TT.
Planck C_l = 2πD_l/(l(l+1)). Model l=1 share of l=1..2500 annotated. Reads/extends A_sky_Cl_l2500.json."""
import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
J = json.load(open("A_sky_Cl_l2500.json")); cl = np.array(J["Cl_model"]); lp = np.array(J["planck_l"]); Dp = np.array(J["planck_Dl_uK2"])
L = np.arange(2, 2501); Cp = 2*np.pi*Dp/(lp*(lp+1))
pm = 100*cl[2:]/cl[2:].sum(); pp = 100*Cp/Cp.sum(); s1 = 100*cl[1]/cl[1:].sum()
cm, cp = np.cumsum(pm), np.cumsum(pp)
reach = lambda c, q: int(L[np.argmax(c >= q)])
R = {f"{q}%": {"model": reach(cm, q), "planck": reach(cp, q)} for q in (50, 90, 99)}
fig, ax = plt.subplots(1, 2, figsize=(14, 5))
ax[0].semilogy(L, pm, color="C3", lw=.8, label="model (blended A/A_all sky)"); ax[0].semilogy(lp, pp, color="k", lw=.8, label="Planck 2018 TT")
ax[0].set_xscale("log"); ax[0].set(xlabel="ℓ", ylabel="% of C_ℓ power over ℓ=2–2500", title="Percent of power per ℓ")
ax[0].text(0.03, 0.04, f"model ℓ=1 (dipole) = {s1:.1f}% of model C_ℓ over ℓ=1–2500", transform=ax[0].transAxes, fontsize=9, bbox=dict(fc="w", ec=".6"))
ax[1].plot(L, cm, color="C3", label="model"); ax[1].plot(lp, cp, color="k", label="Planck")
for q in (50, 90, 99): ax[1].axhline(q, color=".7", lw=.6, ls="--")
ax[1].set_xscale("log"); ax[1].set(xlabel="ℓ", ylabel="cumulative %", title="Cumulative % of C_ℓ power, ℓ=2–2500")
ax[1].text(0.4, 0.1, "\n".join(f"{k}: model ℓ={v['model']}, Planck ℓ={v['planck']}" for k, v in R.items()), transform=ax[1].transAxes, fontsize=9, bbox=dict(fc="w", ec=".6"))
for a in ax: a.grid(alpha=.3); a.legend(fontsize=8)
fig.tight_layout(); fig.savefig("A_sky_Cl_pctpower_vs_planck_l2500.png", dpi=150)
J.update({"model_l1_share_pct": s1, "cumulative_reach_l": R, "pct_model": pm.tolist(), "pct_planck": pp.tolist()}); json.dump(J, open("A_sky_Cl_l2500.json", "w"))
print("l1", s1, R)

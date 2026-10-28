"""ℓ=2..8 power shares (dipole removed): blended A/A_all sky model vs Planck 2018 TT.
Model C_ℓ from A_sky_Cl_blend.json (run A_sky_Cl_blend.py first). Planck: COM_PowerSpect_CMB-TT-full_R3.01.txt,
C_ℓ = 2π D_ℓ/(ℓ(ℓ+1)) — same file/convention as the earlier m-share Planck bars. Writes A_sky_Cl_shares_vs_planck_blend.json/.png."""
import json, numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
J = json.load(open("A_sky_Cl_blend.json")); L = np.arange(2, 9)
d = np.loadtxt("COM_PowerSpect_CMB-TT-full_R3.01.txt"); d = d[(d[:, 0] >= 2) & (d[:, 0] <= 8)]
f = 2*np.pi/(d[:, 0]*(d[:, 0]+1)); Cp = d[:, 1]*f; S = Cp.sum()
pp = 100*Cp/S; lo = 100*d[:, 2]*f/S; hi = 100*d[:, 3]*f/S
out = {"planck_source": "COM_PowerSpect_CMB-TT-full_R3.01.txt", "ell": L.tolist(), "planck_share_pct": pp.tolist(), "model_share_pct": {}}
fig, ax = plt.subplots(figsize=(9, 4.6)); rows = list(J["rows"]); w = 0.8/(len(rows)+1)
for i, k in enumerate(rows):
    c = np.array(J["rows"][k]["Cl"])[2:9]; s = 100*c/c.sum(); out["model_share_pct"][k] = s.tolist()
    ax.bar(L + (i - len(rows)/2)*w, s, w, label="model: " + k)
ax.bar(L + (len(rows)/2)*w, pp, w, color="k", alpha=.75, yerr=[lo, hi], capsize=2, label="Planck 2018 TT")
ax.set(xlabel="ℓ", ylabel="% of C_ℓ power over ℓ=2..8", title="Dipole removed: C_ℓ shares, blended A/A_all sky vs Planck")
ax.legend(fontsize=7); ax.grid(axis="y", alpha=.3); fig.tight_layout(); fig.savefig("A_sky_Cl_shares_vs_planck_blend.png", dpi=150)
json.dump(out, open("A_sky_Cl_shares_vs_planck_blend.json", "w"), indent=1)
print("planck", np.round(pp, 1))
for k, v in out["model_share_pct"].items(): print(k, np.round(v, 1))

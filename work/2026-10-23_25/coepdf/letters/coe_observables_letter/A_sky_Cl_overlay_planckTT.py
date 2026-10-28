"""Blended A/A_all sky (untilted row; rows identical) — monopole removed (healpy.remove_monopole), dipole kept.
D_l overlaid on Planck 2018 TT (COM_PowerSpect_CMB-TT-full_R3.01.txt).
Panel 1: model in uK^2 IF dA/A maps 1:1 to dT/T: D_l * (T_CMB*1e6)^2, T_CMB=2.7255 K.
Panel 2: model scaled so sum_{l=2..30} D_l equals Planck's. Writes A_sky_Cl_overlay_pct.json (shared with pct fig)."""
import json, numpy as np, healpy as hp, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import make_fig05_sky_orientation_values_nobands_blend as V
NS, LMAX, TCMB = 128, 30, 2.7255
m = V.load(V.base.RING_SCRIPTS[0], 0); geo = m.compute_geometry(m.migkas_scan_stops())
n = np.stack(hp.pix2vec(NS, np.arange(hp.nside2npix(NS))), -1)
mp = hp.remove_monopole(V.A_on_dirs(m, geo["D"], n))
cl = hp.anafast(mp, lmax=LMAX); l = np.arange(LMAX + 1); Dm = l*(l+1)*cl/(2*np.pi)
d = np.loadtxt("COM_PowerSpect_CMB-TT-full_R3.01.txt"); d = d[d[:, 0] <= LMAX]
lp, Dp, elo, ehi = d.T
K1 = (TCMB*1e6)**2; Du = Dm*K1
K2 = Dp.sum()/Dm[2:].sum()
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
for a, Dmod, lab, k in [(ax[0], Du, "model, if δA/A ↔ ΔT/T 1:1 (×T_CMB²)", K1), (ax[1], Dm*K2, "model scaled to Planck ℓ=2–30 total", K2)]:
    a.errorbar(lp, Dp, yerr=[elo, ehi], fmt="ko", ms=4, capsize=2, label="Planck 2018 TT")
    a.semilogy(l[1:], Dmod[1:], "o-", color="C3", ms=4, label=lab)
    a.text(0.03, 0.04, f"scale = {k:.4g} μK² per (A/A_all)²", transform=a.transAxes, fontsize=9)
    a.set(xlabel="ℓ"); a.grid(alpha=.3); a.legend(fontsize=8)
ax[0].set(ylabel=r"$D_\ell=\ell(\ell+1)C_\ell/2\pi$  [μK²]", title="Monopole removed, dipole kept: 1:1 temperature mapping")
ax[1].set(title="Model scaled to Planck total power, ℓ=2–30")
fig.tight_layout(); fig.savefig("A_sky_Cl_overlay_planckTT.png", dpi=150)
json.dump({"row": "Untilted, no roll", "nside": NS, "C0_after_remove_monopole": float(cl[0]), "map_mean_after": float(mp.mean()),
           "Cl_model": cl.tolist(), "Dl_model_AA2": Dm.tolist(), "T_CMB_K": TCMB,
           "scale_1to1_uK2": K1, "scale_match_planck_l2_30": K2,
           "Dl_model_uK2_1to1": Du.tolist(), "planck_l": lp.tolist(), "planck_Dl_uK2": Dp.tolist(),
           "planck_err_lo": elo.tolist(), "planck_err_hi": ehi.tolist()}, open("A_sky_Cl_overlay_pct.json", "w"), indent=1)
print("C0", cl[0], "mean", mp.mean(), "K1", K1, "K2", K2)

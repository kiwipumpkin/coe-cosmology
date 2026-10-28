"""Blended A/A_all sky (untilted row), monopole removed, dipole kept; anafast at nside=1024, lmax=2500,
plus nside=512 (lmax=1535) to test resolution dependence. Overlay Planck 2018 TT (COM_PowerSpect_CMB-TT-full_R3.01.txt).
Panel 1: x (T_CMB*1e6)^2, T_CMB=2.7255 K ('if dA/A maps 1:1 to dT/T'); panel 2: scaled so sum D_l(2..2500) = Planck's.
Writes A_sky_Cl_l2500.json."""
import json, numpy as np, healpy as hp, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import make_fig05_sky_orientation_values_nobands_blend as V
TCMB, LMAX = 2.7255, 2500
m = V.load(V.base.RING_SCRIPTS[0], 0); geo = m.compute_geometry(m.migkas_scan_stops())
def spec(ns, lmax):
    n = np.stack(hp.pix2vec(ns, np.arange(hp.nside2npix(ns))), -1)
    mp = hp.remove_monopole(V.A_on_dirs(m, geo["D"], n)); del n
    return hp.anafast(mp, lmax=lmax, use_pixel_weights=False, iter=3)
cl = spec(1024, LMAX); cl512 = spec(512, 1535)
l = np.arange(LMAX + 1); Dm = l*(l+1)*cl/(2*np.pi)
d = np.loadtxt("COM_PowerSpect_CMB-TT-full_R3.01.txt"); d = d[d[:, 0] <= LMAX]; lp, Dp, elo, ehi = d.T
K1 = (TCMB*1e6)**2; K2 = Dp.sum()/Dm[2:].sum()
fig, ax = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
for a, Dmod, lab, k in [(ax[0], Dm*K1, "model, if δA/A ↔ ΔT/T 1:1 (×T_CMB²)", K1), (ax[1], Dm*K2, "model scaled to Planck ℓ=2–2500 total", K2)]:
    a.errorbar(lp, Dp, yerr=[elo, ehi], fmt=".", color="k", ms=2, elinewidth=.4, alpha=.6, label="Planck 2018 TT")
    a.plot(l[1:], Dmod[1:], "-", color="C3", lw=.8, label=lab + " (nside 1024)")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(1, 2600)
    a.text(0.03, 0.04, f"scale = {k:.4g} μK² per (A/A_all)²", transform=a.transAxes, fontsize=9)
    a.set(xlabel="ℓ"); a.grid(alpha=.3); a.legend(fontsize=8, loc="upper right")
ax[0].set(ylabel=r"$D_\ell$ [μK²]", title="Monopole removed, dipole kept: 1:1 temperature mapping")
ax[1].set(title="Model scaled to Planck total D_ℓ power, ℓ=2–2500")
fig.tight_layout(); fig.savefig("A_sky_Cl_overlay_planckTT_l2500.png", dpi=150)
lc = np.array([2, 10, 30, 100, 300, 500, 1000, 1500])
cmp = {int(x): {"Cl_1024": float(cl[x]), "Cl_512": float(cl512[x]), "ratio_512_over_1024": float(cl512[x]/cl[x])} for x in lc}
band = lambda a, b: float(cl512[a:b].sum()/cl[a:b].sum())
json.dump({"row": "Untilted, no roll", "nside": 1024, "lmax": LMAX, "C0": float(cl[0]), "Cl_model": cl.tolist(),
           "Cl_model_nside512": cl512.tolist(), "scale_1to1_uK2": K1, "scale_match_planck_l2_2500": K2,
           "planck_l": lp.tolist(), "planck_Dl_uK2": Dp.tolist(), "planck_err_lo": elo.tolist(), "planck_err_hi": ehi.tolist(),
           "nside_compare": cmp, "band_ratio_512_over_1024": {"2-30": band(2, 31), "31-300": band(31, 301), "301-1000": band(301, 1001), "1001-1535": band(1001, 1536)}},
          open("A_sky_Cl_l2500.json", "w"))
print("C0", cl[0], "K1", K1, "K2", K2); print(json.dumps(cmp, indent=0)); print(json.load(open("A_sky_Cl_l2500.json"))["band_ratio_512_over_1024"])

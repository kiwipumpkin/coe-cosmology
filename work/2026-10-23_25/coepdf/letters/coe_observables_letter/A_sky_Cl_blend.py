"""Spherical power spectrum of the blended A/A_all sky (fig8_orientation_seq_values_nobands_blend, panel B).
Same value formula/orientation as make_fig05_sky_orientation_values_nobands_blend.py (imported):
value(n) = w*A(theta)+(1-w)*A(360-theta), theta=arccos(n.D), w=(1+cos psi)/2, cos psi=(n.y)/sqrt(1-(n.D)^2).
HEALPix nside=128, full sky, no mask, mean removed, healpy.anafast lmax=30. Writes A_sky_Cl_blend.json/.png."""
import json, numpy as np, healpy as hp, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import make_fig05_sky_orientation_values_nobands_blend as V
NS, LMAX = 128, 30
LAB = ["Untilted, no roll", "+15° tilt", "+15° tilt, gradient", "+15° tilt, gradient, +80.8° roll"]
n = np.stack(hp.pix2vec(NS, np.arange(hp.nside2npix(NS))), -1)  # Galactic (cos b cos l, cos b sin l, sin b)
out = {"nside": NS, "lmax": LMAX, "rows": {}}
ell = np.arange(LMAX + 1)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for i, s in enumerate(V.base.RING_SCRIPTS):
    m = V.load(s, i); geo = m.compute_geometry(m.migkas_scan_stops())
    mp = V.A_on_dirs(m, geo["D"], n); mp = mp - mp.mean()
    cl = hp.anafast(mp, lmax=LMAX)
    out["rows"][LAB[i]] = {"script": s, "D": list(map(float, geo["D"])), "Cl": cl.tolist(),
                           "map_rms": float(mp.std())}
    for a, y in zip(ax, [cl, ell*(ell+1)*cl/(2*np.pi)]):
        a.semilogy(ell[1:], y[1:], "o-" if i == 0 else ["s--", "^:", "x-."][i-1], ms=4, label=LAB[i])
ax[0].set(xlabel="ℓ", ylabel=r"$C_\ell$  [(A/A$_{all}$)$^2$]", title="Blended A/A_all sky: C_ℓ (mean removed, full sky)")
ax[1].set(xlabel="ℓ", ylabel=r"$D_\ell=\ell(\ell+1)C_\ell/2\pi$", title="Same, as D_ℓ")
for a in ax: a.legend(fontsize=7); a.grid(alpha=.3)
C = np.array([r["Cl"] for r in out["rows"].values()])
out["max_rel_spread_l1_8"] = float(np.max((C.max(0)-C.min(0))[1:9]/C.mean(0)[1:9]))
fig.tight_layout(); fig.savefig("A_sky_Cl_blend.png", dpi=150)
json.dump(out, open("A_sky_Cl_blend.json", "w"), indent=1)
print("spread", out["max_rel_spread_l1_8"])
for k, r in out["rows"].items(): print(k, ["%.3e" % c for c in r["Cl"][1:9]])

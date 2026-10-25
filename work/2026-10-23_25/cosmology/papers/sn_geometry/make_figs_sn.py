import numpy as np, json, warnings, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from PIL import Image
warnings.filterwarnings("ignore")
import sncore as s
plt.rcParams.update({"font.size": 9, "font.family": "serif"})
def save(fig, name):
    fig.savefig(name + ".png", dpi=200, bbox_inches="tight"); im = Image.open(name + ".png").convert("RGB")
    if im.width > 1700: im = im.resize((1700, int(im.height * 1700 / im.width)), Image.LANCZOS)
    im.save(name + "_web.jpg", quality=88)
S = s.load(0.023); z = S["z"]; COE = s.DIRS["CoE"]
fit = min(json.load(open("fit_z0.023.json")), key=lambda r: r["chi2"])
muL = s.mu_frw(z, 0.327, 0.673); cL, oL = s.chi2_marg(S, muL)
res = S["mb"] - oL - muL; err = np.sqrt(np.diag(S["C"]))
cpsi = S["n"] @ s.unit(*COE)
zg = np.geomspace(0.023, 2.3, 200); muLg = s.mu_frw(zg, 0.327, 0.673)
def curve_frw(Om):
    o = s.chi2_marg(S, s.mu_frw(z, Om, 0.0))[1]; return s.mu_frw(zg, Om, 0.0) + o - muLg - oL
def ltb(D, LD, mc, Om, lc, bc):
    p, t0, ro = s.geom_from(D, LD, mc, Om); T = s.dl_table(p, t0, ro)
    o = s.chi2_marg(S, s.mu_ltb(T, s.dl_nodes(T, z), S["n"], lc, bc))[1]
    n0 = s.unit(lc, bc); DLg = s.dl_nodes(T, zg)
    return {k: s.mu_ltb(T, DLg, np.tile(v, (zg.size, 1)), lc, bc) + o - muLg - oL
            for k, v in [("toward", n0), ("away", -n0), ("perp", np.cross(n0, [0, 0, 1]) / np.linalg.norm(np.cross(n0, [0, 0, 1])))]}
chris = ltb(41.53, 2.919 / 41.53, 0.02395, 0.0, *COE)
bestc = ltb(fit["D"], fit["LD"], fit["mc"], fit["Om"], fit["l"], fit["b"])
edges = np.geomspace(0.023, 2.3, 12)
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.7), sharey=True)
groups = [("toward CoE, $\\cos\\psi>1/3$", cpsi > 1 / 3, "toward"), ("perpendicular, $|\\cos\\psi|<1/3$", np.abs(cpsi) <= 1 / 3, "perp"),
          ("away from CoE, $\\cos\\psi<-1/3$", cpsi < -1 / 3, "away")]
for a, (title, m, key) in zip(ax, groups):
    a.scatter(z[m], res[m], s=1.5, c="0.75", lw=0)
    zc, mu_, e_ = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = m & (z >= lo) & (z < hi)
        if k.sum() < 4: continue
        w = 1 / err[k]**2; zc.append(np.exp(np.mean(np.log(z[k])))); mu_.append(np.sum(w * res[k]) / w.sum()); e_.append(1 / np.sqrt(w.sum()))
    a.errorbar(zc, mu_, e_, fmt="o", ms=3, color="k", lw=0.8, zorder=5)
    a.axhline(0, color="k", lw=0.6, label=r"flat $\Lambda$CDM")
    a.plot(zg, curve_frw(0.0), color="C2", lw=1, label="Milne")
    a.plot(zg, curve_frw(0.05), color="C2", ls="--", lw=1, label=r"open FRW, $\Omega_m=0.05$")
    a.plot(zg, chris[key], color="C3", lw=1.3, label="CoE LTB (Chris geometry)")
    a.plot(zg, bestc["perp"] if key == "perp" else bestc[key], color="C0", lw=0.8, ls=":", label="best free LTB, own axis")
    a.set_xscale("log"); a.set_ylim(-0.35, 0.35); a.set_title(title, fontsize=8); a.set_xlabel("$z$")
ax[0].set_ylabel(r"$\mu-\mu_{\Lambda{\rm CDM}}$ [mag]"); ax[0].legend(fontsize=6, loc="lower left", frameon=False)
fig.tight_layout(); save(fig, "fig_hubble_resid")
# sky map
d = json.load(open("direction.json"))["z0.023"]
fig = plt.figure(figsize=(7.0, 3.8)); a = fig.add_subplot(111, projection="mollweide")
X = lambda l: -np.radians(((np.asarray(l) + 180) % 360) - 180); Y = lambda b: np.radians(b)
lo = z < 0.1
sc = a.scatter(X(S["l"][lo]), Y(S["b"][lo]), c=res[lo], s=5, cmap="coolwarm", vmin=-0.3, vmax=0.3, lw=0)
mk = [("CMB dipole", 264.02, 48.25, "k", "*"), ("Migkas", 280, -20, "m", "s"), ("CoE", 193.17, 0, "r", "D"),
      ("anti-CoE", 13.17, 0, "r", "d"), ("SN dipole $z<0.1$", d["dip_LCDM_0_0.1"]["l"], d["dip_LCDM_0_0.1"]["b"], "b", "P"),
      ("SN dipole all $z$", d["dip_LCDM"]["l"], d["dip_LCDM"]["b"], "c", "P"), ("best free LTB axis", fit["l"], fit["b"], "g", "o")]
for n, l, b, c, m in mk:
    a.plot(X(l), Y(b), m, color=c, ms=8, mec="k", mew=0.6, label=n)
a.set_xticklabels(["150°", "120°", "90°", "60°", "30°", "0°", "330°", "300°", "270°", "240°", "210°"], fontsize=7)
a.grid(alpha=0.3); a.legend(fontsize=6.5, loc="lower center", bbox_to_anchor=(0.5, -0.32), ncol=4, frameon=False)
cb = fig.colorbar(sc, ax=a, shrink=0.6, pad=0.02); cb.set_label(r"$\mu-\mu_{\Lambda{\rm CDM}}$, $z<0.1$ [mag]", fontsize=8)
save(fig, "fig_sky")

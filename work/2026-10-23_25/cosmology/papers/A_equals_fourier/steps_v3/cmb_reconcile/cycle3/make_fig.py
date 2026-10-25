"""Cycle 3 summary figure (no equations in panels)."""
import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import migkas_A as M
from slices import fit_curve, far_af, near_af, SL, TH1
S = json.load(open("slices2.json")); V = json.load(open("validity_ladder.json")); BE = json.load(open("mix3_edge.json"))
NU = json.load(open("mix3_null.json")) + json.load(open("mix3_null2.json")); HS = json.load(open("hotspot.json"))
fig = plt.figure(figsize=(17, 13)); gs = fig.add_gridspec(3, 4, hspace=0.42, wspace=0.3)
names = ["z<0.03", "0.03-0.07", "z>0.07", "all"]
for k, s in enumerate(names):
    ax = fig.add_subplot(gs[0, k]); e = S[s]; A = np.array(e["A"]); er = np.array(e["A_err"]); ok = np.isfinite(A)
    ax.fill_between(M.LC[ok], (A - er)[ok], (A + er)[ok], color="0.85"); ax.plot(M.LC[ok], A[ok], "k.", ms=4, label="Migkas 2020 clusters, rebuilt")
    zm = e["z_median"]; lc, y = M.LC[ok], A[ok]
    _, cf = fit_curve(far_af(299792.458 * zm / 71.6 / 1e3), lc, y); ax.plot(TH1, cf, "C0", label="far core (refit l0+sine) R2 %.2f" % e["far_refit"]["R2"])
    for D, c in ((0.2, "C3"), (0.5, "C2")):
        _, cn = fit_curve(near_af((0.08, D, zm)), lc, y); ax.plot(TH1, cn, c, lw=1, label="near core D=%.1f Gpc R2 %.2f" % (D, e["near"][str(D)]["refit"]["R2"]))
    ax.set_title("%s  (n=%d, median z %.3f)" % (s, e["n"], zm), fontsize=10); ax.set_xlabel("Galactic longitude (deg)"); ax.set_xlim(0, 360)
    if k == 0: ax.set_ylabel("A / A_all")
    ax.legend(fontsize=6.5, loc="lower left")
ax = fig.add_subplot(gs[1, 0]); zz = [S[s]["z_median"] for s in names[:3]]
ax.errorbar(zz, [S[s]["harm"]["amp2"] for s in names[:3]], [S[s]["harm_err"]["amp2"] for s in names[:3]], fmt="ko", capsize=3, label="data: 2l amplitude")
ax.plot(zz, [S[s]["coherent2"] for s in names[:3]], "ks", mfc="none", label="data: part aligned with full-sample phase")
ax.plot(zz, [S[s]["iso_amp2_p95"] for s in names[:3]], "k:", label="isotropic sky, 95th percentile")
for D, c in ((0.2, "C3"), (0.3, "C1"), (0.5, "C2")): ax.plot(zz, [S[s]["near"][str(D)]["model_amp2_nosine"] for s in names[:3]], c + "-o", ms=3, label="near core D=%.1f Gpc (model)" % D)
ax.plot(zz, [S[s]["far_model_amp2_nosine"] for s in names[:3]], "C0-o", ms=3, label="far core headline (model)")
ax.set_xlabel("slice median redshift"); ax.set_ylabel("2l Fourier amplitude of A(l)"); ax.set_title("Is the 2l pattern only nearby?", fontsize=10); ax.legend(fontsize=6.5)
ax = fig.add_subplot(gs[1, 1])
ms = [v["ms"] for v in V]; ax.semilogx(ms, [v["R2_s"] for v in V], "k-", lw=0.8)
for v in V: ax.plot(v["ms"], v["R2_s"], "o", color="C2" if v["ok"] else "C3")
for b in BE:
    if "cb" in b and b["Ab_frac"] > 0.99: ax.plot(b["ms"], b["R2_s"], "^", color="C0", ms=5)
ax.plot([], [], "o", color="C2", label="single bang, regular"); ax.plot([], [], "o", color="C3", label="single bang, NOT regular (1+2E<0 inside)")
ax.plot([], [], "^", color="C0", label="with young-core bang-time bump, regular"); ax.axhline(0.85, color="0.5", ls="--", lw=0.8)
ax.set_xlabel("core mass / headline core mass"); ax.set_ylabel("A(l) R2 at z=0.3"); ax.set_title("Far core, exact GR: valid spacetimes?", fontsize=10); ax.legend(fontsize=6.5)
ax = fig.add_subplot(gs[1, 2])
for w, c in ((10.0, "C0"), (16.5, "C1"), (20.0, "C4"), (24.0, "C5"), (28.0, "C6"), (33.0, "C2"), (66.0, "C3")):
    E = sorted([e for e in NU if e["w"] == w and "cb" in e], key=lambda e: e["Ab_frac"])
    if E: ax.semilogy([-e["Ab_frac"] * 1e5 for e in E], [max(abs(x) for x in e["cb"][2:5]) for e in E], "-o", color=c, ms=3, label="bump width %.1f Gpc" % w)
ax.axhline(1e-5, color="k", ls="--", lw=0.8, label="CMB bound 1e-5")
ax.set_xlabel("older-core bump amplitude (units of 1e-5 of the age)"); ax.set_ylabel("largest of |a2|, |a3|, |a4|")
ax.set_title("Two-mode mix, small core, 48 rays", fontsize=10); ax.legend(fontsize=6)
ax = fig.add_subplot(gs[1, 3])
for b in BE:
    if "cb" in b and b["Ab_frac"] > 0.99: ax.plot(b["R2_s"], max(abs(x) for x in b["cb"][2:5]), "^", color="C0")
    elif "cb" in b: ax.plot(b["R2_s"], max(abs(x) for x in b["cb"][2:5]), "o", color="C2")
for e in NU:
    if "cb" in e: ax.plot(e["R2_s"], max(abs(x) for x in e["cb"][2:5]), ".", color="C1", ms=3)
ax.set_yscale("log"); ax.axhline(1e-5, color="k", ls="--", lw=0.8); ax.axvline(0.85, color="0.5", ls="--", lw=0.8)
ax.set_xlabel("A(l) R2 at z=0.3"); ax.set_ylabel("largest of |a2|, |a3|, |a4|"); ax.set_title("Regular models only: target box is lower right", fontsize=10)
ax = fig.add_subplot(gs[2, 0:2])
for n, c in (("reg_ms0.7_w12", "C0"), ("reg_ms1.0_w13.5", "C1")):
    h = HS[n]; ax.semilogy(h["psi_dense"], np.abs(h["dT_dense"]), color=c, label=n.replace("reg_", "regular: core x").replace("_w", ", bump width ") + " Gpc")
h = HS["headline_single_bang"]; bad = np.array(h["nerr_dense"]) > 1e-6
ax.semilogy(np.array(h["psi_dense"])[~bad], np.abs(np.array(h["dT_dense"]))[~bad], "k.", ms=3, label="single-bang headline (valid rays)")
ax.semilogy(np.array(h["psi_dense"])[bad], np.abs(np.array(h["dT_dense"]))[bad], "rx", ms=3, label="single-bang headline (broken rays)")
ax.axhline(1e-5, color="k", ls="--", lw=0.8); ax.set_xlabel("angle from the CoE direction (deg)"); ax.set_ylabel("|temperature contrast|")
ax.set_title("CoE hot spot, 96 rays over 0.1 to 60 deg", fontsize=10); ax.legend(fontsize=7)
ax = fig.add_subplot(gs[2, 2:4])
for n, c in (("reg_ms0.7_w12", "C0"), ("reg_ms1.0_w13.5", "C1")):
    h = HS[n]; ax.plot(h["psi96"], 1 / (1 + np.array(h["z96"])) * 1090, "o-", ms=2, color=c, label=n)
ax.set_xlabel("angle from the CoE direction (deg), 96 Gauss rays"); ax.set_ylabel("observed temperature / standard (before boost)")
ax.set_yscale("log"); ax.set_title("Whole-sky temperature, regular young-core models", fontsize=10); ax.legend(fontsize=7)
fig.savefig("cycle3_summary.png", dpi=110, bbox_inches="tight"); print("ok")

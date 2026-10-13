import json, sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "../cycle3"); import migkas_A as M
R = json.load(open("solve.json")); S2 = json.load(open("../cycle3/slices2.json")); RT = json.load(open("rotate.json")); TH1 = np.arange(0, 360, 1.0); LC = M.LC
pick = lambda s, D=49.1746, cone=False, mode="median": [r for r in R if r["slice"] == s and abs(r["D"] - D) < 1e-3 and r["mode"] == mode and r["cone"] == cone and r["l0_branch"] == 0][0]
col = {"0.03-0.07": "C1", "z>0.07": "C3", "all": "k"}
# Letter gold is 2x2 (1123x1014): three A/slice-mean panels + degeneracy
fig, axs = plt.subplots(2, 2, figsize=(9.5, 8.6), gridspec_kw=dict(wspace=0.27, hspace=0.32))
axs_flat = axs.ravel()
for ax, s in zip(axs_flat[:3], ("0.03-0.07", "z>0.07", "all")):
    p = np.array(S2[s]["A"], float); e = np.array(S2[s]["A_err"], float); ok = np.isfinite(p); m = np.nanmean(p)
    ax.fill_between(LC[ok], (p - e)[ok] / m, (p + e)[ok] / m, color="0.85"); ax.plot(LC[ok], p[ok] / m, ".", color=col[s], ms=4, label="Migkas 2020 slice (rebuilt)")
    r = pick(s); ax.plot(TH1, r["curve"], color=col[s], lw=2, label="solved model: R2 %.2f" % r["R2_point"])
    rc = pick(s, cone=True); cb = np.convolve(np.r_[np.array(rc["curve"])[-45:], rc["curve"], np.array(rc["curve"])[:45]], np.ones(91) / 91, mode="valid")
    ax.plot(TH1, cb, color=col[s], ls="--", lw=1, label="solved for binned profile, shown binned: R2 %.2f" % rc["R2_conebinned"])
    if s == "0.03-0.07":
        A = np.array(pick("all")["curve"]); dl = RT["all|point"]["dl"]
        ax.plot(TH1, np.interp((TH1 - dl) % 360, TH1, A, period=360), "k:", lw=1.5, label="black model rotated by %.0f deg: R2 %.2f" % (dl, RT["all|point"]["R2"]))
    ax.set_title("%s (median distance %.2f Gpc)" % (s, r["Rs"][0]), fontsize=10); ax.set_xlabel("Galactic longitude (deg)"); ax.set_xlim(0, 360); ax.legend(fontsize=7, loc="lower left")
axs_flat[0].set_ylabel("A / slice mean")
axs_flat[2].set_ylabel("A / slice mean")
ax = axs_flat[3]
for s in ("0.03-0.07", "z>0.07", "all"):
    rr = sorted([r for r in R if r["slice"] == s and r["mode"] == "median" and not r["cone"] and r["l0_branch"] == 0], key=lambda r: r["D"])
    ax.loglog([r["D"] for r in rr], [r["Mc_kg"] for r in rr], "o-", color=col[s], label="%s: core share of critical mass %.4f to %.4f" % (s, min(r["mcc"] for r in rr), max(r["mcc"] for r in rr)))
ax.set_xticks([20,30,49.17,94,200]); ax.set_xticklabels(["20","30","49","94","200"]); ax.minorticks_off(); ax.set_xlabel("CoE distance D (Gpc), not constrained by the data"); ax.set_ylabel("solved core mass (kg)"); ax.set_title("Degeneracy: every D fits exactly", fontsize=10); ax.legend(fontsize=7)
fig.savefig("termslice_summary.png", dpi=140, bbox_inches="tight"); print("ok")

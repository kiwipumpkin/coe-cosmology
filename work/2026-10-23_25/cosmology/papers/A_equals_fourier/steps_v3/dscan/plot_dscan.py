import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
N = json.load(open("dscan_newt.json")); X = json.load(open("dscan_newt_extra.json")); Rl = json.load(open("dscan_rel.json")); C = json.load(open("dscan_rel_checks.json"))
def ser(rows, key="fit", val=None):
    rows = sorted([r for r in rows if val is None or r[key] == val], key=lambda r: r["DR"]); return np.array([r["DR"] for r in rows]), np.array([r["R2"] for r in rows]), np.array([100 * r["share"] for r in rows])
allN = N + X["extra"]
fr = ser(allN, val="free"); om = ser(allN, val="omb")
thin = [r for r in Rl if r["kout"] < 2] + [dict(DR=10, R2=C[1]["R2"], share=C[1]["share"])]
thick = [r for r in Rl if r["kout"] > 2] + [dict(DR=20, R2=C[6]["R2"], share=C[6]["share"])]
lt = ser(thin); lk = ser(thick)
os49 = ser([r for L in X["os_sens"] for r in L if r["Os"] == 0.049])
fig, ax = plt.subplots(1, 2, figsize=(13, 6), dpi=110)
st = [("Newtonian, shell density free (best fit at each D)", fr, "C0", "o-"), ("Newtonian, shell density free, slice at 0.049", os49, "C0", "s:"),
      ("Newtonian, shell density 0.049", om, "C1", "o-"), ("Relativistic LTB, thin-shell branch", lt, "C3", "o-"), ("Relativistic LTB, thick-shell branch", lk, "C3", "D")]
for lab, s, col, m in st:
    ax[0].plot(s[0], s[1], m, color=col, label=lab, ms=6, mfc="white" if m == "D" else None)
    ax[1].plot(s[0], s[2], m, color=col, label=lab, ms=6, mfc="white" if m == "D" else None)
for h, col in ((0.909, "C0"), (0.893, "C1"), (0.885, "C3")): ax[0].axhspan(h - 0.01, h, color=col, alpha=0.08)
for y in (10, 50): ax[1].axhline(y, color="gray", ls="--", lw=0.8); ax[1].text(1.9, y * 1.08, f"{y}%", color="gray")
for a in ax: a.set_xscale("log"); a.set_xlabel("D / R   (R = 14.3 Gpc held fixed)"); a.grid(alpha=0.3, which="both")
ax[0].set_ylabel("R² vs Migkas A/A_all (with sine term)"); ax[0].set_ylim(0.78, 0.92); ax[0].set_title("Fit quality against D (shaded: within 0.01 of each headline)")
ax[1].set_yscale("log"); ax[1].set_ylabel("Core share of total mass (%)"); ax[1].set_title("Core share against D"); ax[1].set_ylim(0.02, 150)
h, l = ax[1].get_legend_handles_labels(); fig.legend(h, l, fontsize=9, loc="lower center", ncol=3, frameon=False)
sec = ax[0].secondary_xaxis("top", functions=(lambda x: x * 14.3, lambda x: x / 14.3)); sec.set_xlabel("D (Gpc)")
sec = ax[1].secondary_xaxis("top", functions=(lambda x: x * 14.3, lambda x: x / 14.3)); sec.set_xlabel("D (Gpc)")
fig.tight_layout(rect=(0, 0.1, 1, 1)); fig.savefig("dscan_R2_coreshare.jpg", pil_kwargs=dict(quality=88))

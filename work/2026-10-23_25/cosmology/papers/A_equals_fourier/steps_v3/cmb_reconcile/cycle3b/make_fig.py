import json, sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "../cycle3"); import migkas_A as M
O = json.load(open("compare.json")); S2 = json.load(open("../cycle3/slices2.json")); TH1 = np.arange(0, 360, 1.0); LC = M.LC
names = ["z<0.03", "0.03-0.07", "z>0.07", "all"]; col = {"free": "C0", "omb": "C1"}
lab = {"free": "Newtonian, shell density free (D 94 Gpc)", "omb": "Newtonian, shell density 0.049 (D 49 Gpc)"}
fig = plt.figure(figsize=(17, 10)); gs = fig.add_gridspec(2, 4, hspace=0.35, wspace=0.28)
for k, s in enumerate(names):
    ax = fig.add_subplot(gs[0, k]); d = np.array(S2[s]["A"]); er = np.array(S2[s]["A_err"]); ok = np.isfinite(d)
    ax.fill_between(LC[ok], (d - er)[ok], (d + er)[ok], color="0.85"); ax.plot(LC[ok], d[ok], "k.", ms=4, label="Migkas 2020 slice, rebuilt")
    for n in ("free", "omb"):
        e = O[n][s]; ax.plot(TH1, e["curve_point"], col[n], label="%s: R2 %.2f" % (n, e["R2_point"]))
        ax.plot(LC, e["curve_percluster"], col[n], ls="--", lw=0.8, label="%s, per-cluster binned: R2 %.2f" % (n, e["R2_percluster"]))
    ax.set_title("%s (n=%d, median distance %.2f Gpc)" % (s, O["omb"][s]["n"], O["omb"][s]["r_median_Gpc"]), fontsize=9.5)
    ax.set_xlabel("Galactic longitude (deg)"); ax.set_xlim(0, 360); ax.legend(fontsize=6.3, loc="lower left")
    if k == 0: ax.set_ylabel("A / A_all")
ax = fig.add_subplot(gs[1, 0:2])
for n in ("free", "omb"):
    C = O["continuous"][n]; R = [c["Rs"] for c in C]
    ax.semilogx(R, [c["max"] for c in C], col[n], label="%s: predicted maximum" % n); ax.semilogx(R, [c["min"] for c in C], col[n], ls="--", label="%s: predicted minimum" % n)
for s, mk in zip(names, "osD^"):
    e = O["omb"][s]["data"]; r = O["omb"][s]["r_median_Gpc"]; lo, hi = O["omb"][s]["r_range_Gpc"]
    ax.errorbar(r, e["max"], xerr=[[r - lo], [hi - r]], fmt="k" + mk, ms=7, capsize=2, label="data %s: max (filled) / min (open)" % s)
    ax.plot(r, e["min"], "k" + mk, mfc="none", ms=7)
ax.axvline(14.3, color="0.6", lw=0.8, ls=":"); ax.text(14.3, 5, " headline sampling sphere", fontsize=7, color="0.4", rotation=90)
ax.set_ylim(0, 470); ax.set_yticks(range(0, 361, 60)); ax.set_xlabel("sampling distance from Earth (Gpc)"); ax.set_ylabel("longitude of A(l) max / min (deg)")
ax.set_title("Predicted extremum positions vs sampling distance (full-sample fits held fixed)", fontsize=10); ax.legend(fontsize=6.3, ncol=3, loc="upper center")
ax = fig.add_subplot(gs[1, 2])
for n in ("free", "omb"):
    C = O["continuous"][n]; ax.semilogx([c["Rs"] for c in C], [c["ph2"] for c in C], col[n], label="%s: predicted 2l phase" % n)
for s, mk in zip(names, "osD^"):
    h = S2[s]["harm"]; he = S2[s]["harm_err"]; ph = O["omb"][s]["data"]["phase2"]; perr = np.degrees(0.5 * he["amp2"] / h["amp2"])
    ax.errorbar(O["omb"][s]["r_median_Gpc"], ph, yerr=perr, fmt="k" + mk, ms=7, capsize=3, label="data %s" % s)
ax.set_ylim(0, 180); ax.set_xlabel("sampling distance from Earth (Gpc)"); ax.set_ylabel("2l phase (deg, modulo 180)")
ax.set_title("2l phase vs sampling distance", fontsize=10); ax.legend(fontsize=6.5)
ax = fig.add_subplot(gs[1, 3]); x = np.arange(4)
for j, n in enumerate(("free", "omb")):
    ax.plot(x + 0.08 * (j - 0.5), [O[n][s]["R2_point"] for s in names], "o", color=col[n], label="%s, point curve" % n)
    ax.plot(x + 0.08 * (j - 0.5), [O[n][s]["R2_percluster"] for s in names], "x", color=col[n], label="%s, per-cluster binned" % n)
ax.axhline(0, color="k", lw=0.6); ax.axhline(0.85, color="0.6", ls="--", lw=0.8); ax.set_xticks(x); ax.set_xticklabels(names, fontsize=8)
ax.set_ylabel("R2 (no refitting)"); ax.set_title("R2 per slice at the slice's own distances", fontsize=10); ax.legend(fontsize=6.5)
fig.savefig("cycle3b_summary.png", dpi=110, bbox_inches="tight"); print("ok")

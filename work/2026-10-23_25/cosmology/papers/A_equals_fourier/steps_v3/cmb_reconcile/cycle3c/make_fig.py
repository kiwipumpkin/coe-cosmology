import json, sys, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "../cycle3"); import migkas_A as M
R = json.load(open("noise_test.json")); LC = M.LC; names = ["z<0.03", "0.03-0.07", "z>0.07", "all"]
fig, axs = plt.subplots(2, 4, figsize=(18, 10), gridspec_kw=dict(hspace=0.75, wspace=0.25))
for k, s in enumerate(names):
    e = R[s]; d = np.array(e["data"], float); lo, hi = np.array(e["boot_band"][0], float), np.array(e["boot_band"][1], float); ok = np.isfinite(d)
    ax = axs[0, k]; ax.fill_between(LC[ok], lo[ok], hi[ok], color="0.82", label="data bootstrap band (16-84%)"); ax.plot(LC[ok], d[ok], "k.", ms=4, label="data (rebuilt Migkas 2020)")
    for n, c in (("omb", "C1"), ("free", "C0")):
        m = e["models"][n]; mm = np.array(m["model"], float); mm[~ok] = np.nan
        ax.plot(LC, mm, c, lw=1.8, label="%s at slice distances: chi2 %.1f for 6, p %.2g (mock-cal. %.2g)" % (n, m["chi"]["chi2_F"], m["chi"]["p_F"], m["p_F_mockcal"]))
    fi = e["models"]["iso"]; ax.axhline(1, color="0.4", ls=":", lw=1, label="isotropic: chi2 %.1f for 6, p %.2g (mock-cal. %.2g)" % (fi["chi"]["chi2_F"], fi["chi"]["p_F"], fi["p_F_mockcal"]))
    ax.set_title("%s  (n=%d, about %.1f independent points)" % (s, e["n"], e["neff"]), fontsize=10); ax.set_xlim(0, 360); ax.set_xlabel("Galactic longitude (deg)")
    if k == 0: ax.set_ylabel("A / slice mean")
    ax.legend(fontsize=6.4, loc="upper center", bbox_to_anchor=(0.5, -0.16))
    ax = axs[1, k]; bins = np.arange(0, 181, 7.5)
    for n, c in (("omb", "C1"), ("iso", "0.5")):
        ax.hist(np.array(e["models"][n]["mock_ph2"]) % 180, bins=bins, histtype="step", color=c, lw=1.5, density=True, label="mocks: %s" % ("model omb + real noise" if n == "omb" else "isotropic + real noise"))
    ax.axvline(e["data_feats"]["ph2"], color="k", lw=2, label="observed %.0f deg" % e["data_feats"]["ph2"])
    ax.axvline(e["models"]["omb"]["pred"]["ph2"], color="C1", ls="--", label="predicted %.0f deg" % e["models"]["omb"]["pred"]["ph2"])
    mo = e["models"]["omb"]; ax.set_title("2l phase: chance of a miss this large %.2f" % mo["p_ph2"], fontsize=10)
    ax.set_xlabel("2l phase (deg, modulo 180)"); ax.set_xlim(0, 180); ax.legend(fontsize=6.5)
    if k == 0: ax.set_ylabel("mock density")
fig.savefig("cycle3c_summary.png", dpi=110, bbox_inches="tight"); print("ok")

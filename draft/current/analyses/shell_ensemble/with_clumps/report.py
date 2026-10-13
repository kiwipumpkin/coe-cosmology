import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":15})
c=json.load(open("counts.json")); bj=json.load(open("baselines.json")); C=c["counts"]
LBL={"a":"(a) 1e48–1e49","b":"(b) 3e49–3e50","c":"(c) 1e51–1e52"}
unc={0:bj["base0"]}; 
for k in (1,2,3,4): unc[k]=bj["baselines"][str(k)]["best"]["rms_uncomp"]
names=["","2 + 3 (pair at ℓ 75°/315°, just beyond D+d)"]
with open("summary_table.md","w") as f:
    f.write("# Random surround+gradient web on top of k placed clumps (bare, 1,000 layouts per row)\n\n")
    f.write("Clumps: the fitted four (scan.json sym N=4; 7.52e51, 7.52e51, 7.79e51, 7.79e51 kg; all 24.5 Gpc from us in the ring plane). "
            "Clumps 0,1 (ℓ 162°/228°) sit 39.9 Gpc from the CoE (inside D−d); clumps 2,3 (ℓ 75°/315°) sit 73.6 Gpc from the CoE (just beyond D+d).\n\n")
    f.write("Compensation (adopted): each placed clump's mass is removed from the surround region (14.3–42.9 Gpc from us, CoE distance ≥ D−d) in proportion to the exp-lock density ρ0·e^(−r/L), "
            "so the smooth profile is scaled by (1−ε) and keeps its shape/gradient; net added mass = 0. ε ≈ 0.1% for all four.\n\n")
    f.write("|k|clumps placed (best subset of the fitted four)|baseline, uncompensated|baseline, compensated (used)|\n|-|-|-|-|\n")
    for k in range(5):
        b=c["baselines"][str(k)]; f.write(f"|{k}|{b['subset'] if b['subset'] else 'none'}|{unc[k]:.1f}×|{b['rms_over_floor']:.1f}×|\n")
    f.write("\n## Helped / hurt (out of 1,000) relative to each baseline\n\n|N|level|k=0|k=1|k=2|k=3|k=4|\n|-|-|-|-|-|-|-|\n")
    for N in (200,500):
        for lev in "abc":
            f.write(f"|{N}|{LBL[lev]}|"+"|".join(f"{C[f'k{k}_N{N}_{lev}']['placed']['rms_lower']} / {C[f'k{k}_N{N}_{lev}']['placed']['rms_higher']}" for k in range(5))+"|\n")
    for k in range(5):
        base=c["baselines"][str(k)]["rms_over_floor"]
        f.write(f"\n### k = {k} (baseline {base:.1f}×)\n\n|N|level|median RMS/floor|5th pct|min|best-scaled median|best 5th pct|best min|RMS lower / higher (placed)|m3–5 lower|m1–2 higher|3–5 lower & 1–2 not higher|each of m3,m4,m5 lower & each of m1,m2 not higher|best: 3–5 lower|best: 1–2 higher|best: both|\n|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|\n")
        for N in (200,500):
            for lev in "abc":
                X=C[f"k{k}_N{N}_{lev}"]; p=X["placed"]; b=X["best"]
                f.write(f"|{N}|{LBL[lev]}|{p['median']:.1f}|{p['p5']:.1f}|{p['min']:.1f}|{b['median']:.1f}|{b['p5']:.1f}|{b['min']:.1f}|{p['rms_lower']} / {p['rms_higher']}|{p['m35_lower']}|{p['m12_higher']}|{p['m35_lower_and_m12_not_higher']}|{p['each_m3_m4_m5_lower_and_each_m1_m2_not_higher']}|{b['m35_lower']}|{b['m12_higher']}|{b['m35_lower_and_m12_not_higher']}|\n")
print(open("summary_table.md").read())

# Figure 1: helped vs hurt, level b, N=200
ks=range(5); h=[C[f"k{k}_N200_b"]["placed"]["rms_lower"] for k in ks]; base=[c["baselines"][str(k)]["rms_over_floor"] for k in ks]
x=np.arange(5); w=0.38
fig,ax=plt.subplots(figsize=(13,7.5))
ax.bar(x-w/2,h,w,color="#2a9d4a",label="helped: lowered the leftover"); ax.bar(x+w/2,[1000-v for v in h],w,color="#d1495b",label="hurt: raised the leftover")
for xi,v in zip(x,h):
    ax.text(xi-w/2,v+12,str(v),ha="center",fontsize=14,color="#1d6e34",weight="bold"); ax.text(xi+w/2,1000-v+12,str(1000-v),ha="center",fontsize=14,color="#9b2335",weight="bold")
ax.set_xticks(x); ax.set_xticklabels([f"{k} clump{'s' if k!=1 else ''}\n(start {b:.0f}×)" for k,b in zip(ks,base)])
ax.axhline(500,color="0.5",ls=":",lw=1.2); ax.set_ylim(0,1080)
ax.set_ylabel("random layouts (out of 1,000)"); ax.set_xlabel("number of fitted clumps placed first (and the leftover they leave, × the 10⁻⁵ target)")
ax.set_title("After placing 0–4 of the fitted clumps, does a random web of 200 clusters + voids\n(3×10⁴⁹–3×10⁵⁰ kg each, all around us) lower or raise what is left?",fontsize=16)
ax.legend(loc="upper left",fontsize=13); ax.grid(axis="y",alpha=.3)
fig.tight_layout(); fig.savefig("helped_vs_hurt_with_clumps.png",dpi=110); plt.close(fig)

# Figure 2: baseline vs k with random-layout median and 5-95% range (level b, N=200)
med=[C[f"k{k}_N200_b"]["placed"]["median"] for k in ks]; p5=[C[f"k{k}_N200_b"]["placed"]["p5"] for k in ks]; p95=[C[f"k{k}_N200_b"]["placed"]["p95"] for k in ks]
mn=[C[f"k{k}_N200_b"]["placed"]["min"] for k in ks]
fig,ax=plt.subplots(figsize=(12,7.5))
ax.fill_between(x,p5,p95,color="0.8",label="random web added: middle 90% of layouts")
ax.plot(x,med,"o-",color="0.35",lw=2,ms=8,label="random web added: typical (median) layout")
ax.plot(x,mn,"v",color="#2a9d4a",ms=11,label="random web added: best single layout")
ax.plot(x,base,"s-",color="k",lw=3,ms=11,label="clumps only (no random web)")
for xi,b in zip(x,base): ax.text(xi+0.08,b-12,f"{b:.0f}×",fontsize=14,weight="bold")
ax.set_xticks(x); ax.set_xticklabels([f"{k}" for k in ks]); ax.set_xlabel("number of fitted clumps placed")
ax.set_ylabel("leftover CMB pattern (× the 10⁻⁵ target)"); ax.set_ylim(0,max(p95)*1.08); ax.grid(alpha=.3)
ax.set_title("Leftover after placing 0–4 clumps, and what a random web of 200 clusters + voids\n(3×10⁴⁹–3×10⁵⁰ kg each) does on top",fontsize=16)
ax.legend(loc="lower left",fontsize=12.5); fig.tight_layout(); fig.savefig("baseline_vs_k.png",dpi=110); plt.close(fig)
print("figs ok")

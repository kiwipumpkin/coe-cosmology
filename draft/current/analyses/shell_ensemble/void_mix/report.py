import json, glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":15})
O={}
for f in sorted(glob.glob("part_k*.json")): O.update(json.load(open(f)))
json.dump(dict(method="void deficits rescaled so total void deficit = total cluster mass (zero net contrast vs exp-lock profile); "
               "void depth |delta| = 0.5 x rescale factor, capped at 1 (empty); past the cap voids grow in radius; f_v random: uniform 0-1 per layout, "
               "at least 1 void and 1 cluster",counts=O),open("counts.json","w"),indent=1)
FV=["10","25","50","75","90","rand"]; FVL={"rand":"random 0–100%"}; LBL={"a":"(a) 1e48–1e49","b":"(b) 3e49–3e50","c":"(c) 1e51–1e52"}
with open("summary_table.md","w") as f:
    f.write("# Void share of the random web varied (surround + exp-lock gradient, compensated clumps, bare, 1,000 layouts per row)\n\n")
    f.write("Zero net contrast kept: total void deficit = total cluster mass in every layout. Fewer voids → each void carries more deficit: it deepens from δ = −0.5 "
            "down to δ = −1 (completely empty, the cap), and beyond that it gets larger. More voids → each is shallower. "
            "'random' = each layout draws its void share uniformly from 0–100% (at least 1 void and 1 cluster). f_v = 50% reproduces the with_clumps run.\n\n")
    f.write("## Helped / hurt out of 1,000 (relative to the k-clump baseline)\n\n")
    for N in (200,500):
        for lev in "abc":
            f.write(f"\n### N = {N}, level {LBL[lev]}\n\n|k (baseline)|"+"|".join(f"voids {FVL.get(v,v+'%')}" for v in FV)+"|\n|-|"+"-|"*len(FV)+"\n")
            for k in range(5):
                x=O.get(f"k{k}_N{N}_{lev}_fv10")
                if not x: continue
                f.write(f"|{k} ({x['baseline']:.1f}×)|"+"|".join(f"{O[f'k{k}_N{N}_{lev}_fv{v}']['placed']['rms_lower']} / {O[f'k{k}_N{N}_{lev}_fv{v}']['placed']['rms_higher']}" for v in FV)+"|\n")
    f.write("\n## Full detail\n\n|k|N|level|void share|baseline|helped / hurt|median|5th pct|best as drawn|best-scaled median|best-scaled min|3–5 lower & 1–2 not higher|each m3–5 lower & each m1–2 not higher|void depth (median |δ|)|share of layouts at δ=−1 cap|largest void radius (Gpc)|\n|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|\n")
    for k in range(5):
        for N in (200,500):
            for lev in "abc":
                for v in FV:
                    x=O.get(f"k{k}_N{N}_{lev}_fv{v}")
                    if not x: continue
                    p=x["placed"]; b=x["best"]
                    f.write(f"|{k}|{N}|{lev}|{FVL.get(v,v+'%')}|{x['baseline']:.1f}|{p['rms_lower']} / {p['rms_higher']}|{p['median']:.1f}|{p['p5']:.1f}|{p['min']:.1f}|{b['median']:.1f}|{b['min']:.1f}|"
                            f"{p['m35_lower_and_m12_not_higher']}|{p['each_m3_m4_m5_lower_and_each_m1_m2_not_higher']}|{x['void_depth_median']:.2f}|{x['void_depth_capped_frac']:.2f}|{x['max_void_radius_gpc']:.2f}|\n")
print(open("summary_table.md").read()[:6000])

fig,axs=plt.subplots(1,2,figsize=(17,7.5),sharey=True)
x=np.arange(len(FV)); w=0.38
for a,k,t in zip(axs,(0,4),("No clumps placed (start 237×)","All four clumps placed (start 89×)")):
    h=[O[f"k{k}_N200_b_fv{v}"]["placed"]["rms_lower"] for v in FV]
    a.bar(x-w/2,h,w,color="#2a9d4a",label="helped: lowered the leftover"); a.bar(x+w/2,[1000-v for v in h],w,color="#d1495b",label="hurt: raised the leftover")
    for xi,v in zip(x,h):
        a.text(xi-w/2,v+12,str(v),ha="center",fontsize=13,color="#1d6e34",weight="bold"); a.text(xi+w/2,1000-v+12,str(1000-v),ha="center",fontsize=13,color="#9b2335",weight="bold")
    a.set_xticks(x); a.set_xticklabels([("random" if v=="rand" else v+"%") for v in FV]); a.set_title(t,fontsize=16)
    a.set_xlabel("share of the random objects that are voids"); a.axhline(500,color="0.5",ls=":",lw=1.2); a.grid(axis="y",alpha=.3); a.set_ylim(0,1080)
axs[0].set_ylabel("random layouts (out of 1,000)"); axs[0].legend(loc="upper left",fontsize=12.5)
fig.suptitle("Does changing the void / cluster mix help? 200 random objects, 3×10⁴⁹–3×10⁵⁰ kg each, net mass zero",fontsize=17)
fig.tight_layout(rect=(0,0,1,.94)); fig.savefig("helped_vs_hurt_void_mix.png",dpi=110); plt.close(fig); print("fig ok")

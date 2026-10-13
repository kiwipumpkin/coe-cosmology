import json, glob, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False})
D="#d62728"; S="#1f77b4"
sw=json.load(open("swap.json")); nf=json.load(open("near_fit.json")); old=json.load(open("../results.json"))
near=[r["left"] for r in nf["runs"]]; best=min(near)
# 1 near fix
lab=["start","four clumps","deep voids,\nfitted","deep layout,\nvoids made shallow\n(no refit)","shallow voids,\nrefitted"]
val=[237.26,88.81,old["frozen"],sw["shallow_voids_no_clusters"],best]; col=["#9e9e9e","#9e9e9e",D,"#ff9896",S]
fig,ax=plt.subplots(figsize=(8.5,4.6)); ax.bar(range(5),val,color=col,width=.6)
for i,v in enumerate(val): ax.text(i,v+4,f"{v:.0f}×",ha="center")
ax.set_xticks(range(5)); ax.set_xticklabels(lab,fontsize=9); ax.set_ylabel("leftover (× floor)"); ax.set_ylim(0,270)
ax.set_title("Near-region fix (14.3–30 Gpc): deep vs shallow voids"); fig.tight_layout(); fig.savefig("near_fix_deep_vs_shallow.png",dpi=140); plt.close()
# 2 outer change
o=old["outer"]; fz=old["frozen"]
dd=[o["foam-like (grown)"]["median"]-fz,o["random packing"]["median"]-fz]
so=[json.load(open(f"outer_{c}.json")) for c in ("foam","random")]; ds=[np.median([r["left"] for r in s["rows"]])-s["frozen"] for s in so]
fig,ax=plt.subplots(figsize=(7,4.2)); x=np.arange(2)
ax.bar(x-.18,dd,.36,color=D,label=f"deep voids (near frozen at {fz:.1f}×)"); ax.bar(x+.18,ds,.36,color=S,label=f"shallow voids (near frozen at {so[0]['frozen']:.1f}×)")
for i in range(2): ax.text(i-.18,dd[i]+.04,f"+{dd[i]:.1f}×",ha="center"); ax.text(i+.18,ds[i]+.04,f"+{ds[i]:.1f}×",ha="center")
ax.set_xticks(x); ax.set_xticklabels(["foam-like fill","random fill"]); ax.set_ylabel("change in leftover from filling 30–42.9 Gpc"); ax.set_ylim(0,2.5)
ax.legend(frameon=False,fontsize=9); ax.set_title("How much the outer shell moves our view"); fig.tight_layout(); fig.savefig("outer_deep_vs_shallow.png",dpi=140); plt.close()
# 3 growth
B=old["B"]; rules=[("foam","foam_rim_Rv1.5","foam"),("random","random_Rv1.5","random"),("CMB-aware","cmb_aware_Rv1.5","cmb_aware")]
g={r:[json.load(open(f)) for f in sorted(glob.glob(f"grow_{r}_s*.json"))] for _,_,r in rules}
fig,axs=plt.subplots(1,2,figsize=(11,4.5))
for ax,key,ttl in ((axs[0],"min","lowest point while growing"),(axs[1],"final","when the shell is full")):
    dv=[np.median(B[k][key]) for _,k,_ in rules]; sv=[np.median([d[key] for d in g[r]]) for _,_,r in rules]
    x=np.arange(3); ax.bar(x-.18,dv,.36,color=D,label="deep voids + compact clusters"); ax.bar(x+.18,sv,.36,color=S,label="shallow voids, smooth compensation")
    for i in range(3): ax.text(i-.18,dv[i]+15,f"{dv[i]:.0f}×",ha="center",fontsize=9); ax.text(i+.18,sv[i]+15,f"{sv[i]:.0f}×",ha="center",fontsize=9)
    ax.axhline(237.26,color="gray",ls="--"); ax.axhline(88.81,color="#2ca02c",ls="--")
    ax.set_ylim(0,max(dv+sv)*1.15); ax.set_xticks(x); ax.set_xticklabels([a for a,_,_ in rules]); ax.set_title(ttl); ax.set_ylabel("leftover (× floor)")
h,l=axs[0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",ncol=2,frameon=False,fontsize=9); axs[1].text(2.5,95,"four clumps 89×",ha="right",color="#2ca02c",fontsize=9)
fig.suptitle("Growing the shell one void at a time: deep vs shallow (dashed grey = start 237×)"); fig.tight_layout(rect=(0,0.07,1,1)); fig.savefig("growth_deep_vs_shallow.png",dpi=140); plt.close()
json.dump(dict(near_runs=near,best_near=best,swap=sw,outer_shallow={c:dict(frozen=s["frozen"],median=float(np.median([r["left"] for r in s["rows"]])),units=[r["units"] for r in s["rows"]],min_gap=min(r["min_gap"] for r in s["rows"])) for c,s in zip(("foam","random"),so)},
     outer_deep_change=dict(foam=dd[0],random=dd[1]),
     growth_shallow={r:[{k:v for k,v in d.items() if k!="traj"} for d in g[r]] for _,_,r in rules},
     growth_deep={k:dict(min=B[k]["min"],final=B[k]["final"]) for _,k,_ in rules}),open("results_shallow.json","w"),indent=1)
print("ok", [r["n"] for r in nf["runs"]], nf["best"]["amps"])

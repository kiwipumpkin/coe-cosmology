import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":15}); F=1e-5
VM=json.load(open("../counts.json"))["counts"]; FV=["10","25","50","75","90","rand"]; out={}
# ---------- 1) helped by k (line chart)
fig,ax=plt.subplots(figsize=(12,7.5)); x=np.arange(len(FV)); cols=["0.25","#1f77b4","#ff7f0e","#9467bd","#2a9d4a"]
out["helped_by_k_N200_b"]={}
ends={k:VM[f"k{k}_N200_b_fv{FV[-1]}"]["placed"]["rms_lower"] for k in range(5)}
order=sorted(ends,key=lambda k:-ends[k]); ylab={}; prev=None
for k in order:
    y=ends[k] if prev is None else min(ends[k],prev-20); ylab[k]=y; prev=y
for k in range(5):
    h=[VM[f"k{k}_N200_b_fv{v}"]["placed"]["rms_lower"] for v in FV]; base=VM[f"k{k}_N200_b_fv10"]["baseline"]
    out["helped_by_k_N200_b"][k]=dict(baseline=base,helped=dict(zip(FV,h)))
    ax.plot(x,h,"o-",lw=2.5,ms=8,color=cols[k]); ax.text(x[-1]+0.15,ylab[k],("no clumps" if k==0 else f"{k} clump{'s' if k!=1 else ''}")+f" (start {base:.0f}×)",va="center",fontsize=13,color=cols[k],weight="bold")
ax.set_xticks(x); ax.set_xticklabels([("random" if v=="rand" else v+"%") for v in FV]); ax.set_xlim(-0.3,7.6)
ax.axhline(500,color="0.5",ls=":",lw=1.2); ax.text(-0.25,510,"half of the layouts (500)",fontsize=11,color="0.4")
ax.set_ylim(0,560); ax.set_xlabel("share of the random objects that are voids"); ax.set_ylabel("random layouts that helped (out of 1,000)")
ax.set_title("How many of 1,000 random layouts lowered the leftover, for 0–4 placed clumps\n(200 random objects, 3×10⁴⁹–3×10⁵⁰ kg each)",fontsize=16)
ax.grid(alpha=.3); fig.tight_layout(); fig.savefig("helped_vs_hurt_by_k.png",dpi=110); plt.close(fig)

# ---------- 2) masses only vs masses+voids
FM=json.load(open("fit_masses_voids.json")); EN=json.load(open("ensemble_fitted.json")); base0=VM["k0_N200_b_fv50"]["baseline"]
mo=[base0]+[FM[f"k{k}_masses_only"]["rms_over_floor"] for k in (1,2,3,4)]; mv=[base0]+[FM[f"k{k}_masses_voids"]["rms_over_floor"] for k in (1,2,3,4)]
out["fits"]={k:dict(masses_only=mo[k],masses_voids=mv[k]) for k in range(5)}
fig,ax=plt.subplots(figsize=(12,7.5)); kx=np.arange(5)
ax.plot(kx,mo,"s-",color="#d1495b",lw=3,ms=11,label="masses only (clusters)"); ax.plot(kx,mv,"o-",color="#1f77b4",lw=3,ms=11,label="masses or voids (each object free to be either)")
for k in kx:
    ax.text(k+0.07,mo[k]+6,f"{mo[k]:.0f}×",color="#d1495b",fontsize=13,weight="bold")
    if abs(mv[k]-mo[k])>1: ax.text(k+0.07,mv[k]-14,f"{mv[k]:.0f}×",color="#1f77b4",fontsize=13,weight="bold")
ax.axhline(83.4,color="0.5",ls="--",lw=1.2); ax.text(0.02,87,"earlier four-clump fit (83×)",fontsize=11,color="0.4")
ax.set_xticks(kx); ax.set_xlabel("number of objects placed on purpose"); ax.set_ylabel("leftover CMB pattern (× the 10⁻⁵ target)")
ax.set_ylim(0,260); ax.grid(alpha=.3); ax.legend(loc="upper right",fontsize=13)
ax.set_title("Leftover after placing 1–4 objects on purpose: masses only vs masses or voids",fontsize=16)
fig.tight_layout(); fig.savefig("baseline_masses_vs_masses_voids.png",dpi=110); plt.close(fig)

# ---------- 3) how much the helpful layouts helped (k=0, N=200, b)
def stats(fvs):
    r=np.concatenate([np.load(f"../raw_k0_N200_b_fv{v}.npz")["rms"] for v in fvs])/F
    rb=np.concatenate([np.load(f"../raw_k0_N200_b_fv{v}.npz")["rms_best"] for v in fvs])/F
    hel=r<base0; red=100*(base0-r[hel])/base0; redb=100*(base0-rb[hel])/base0; redb_all=100*(base0-rb)/base0
    P=lambda a:dict(median=float(np.median(a)),p75=float(np.percentile(a,75)),p90=float(np.percentile(a,90)),best=float(a.max()))
    return dict(n_layouts=len(r),n_helped=int(hel.sum()),as_drawn_pct_reduction=P(red),as_drawn_leftover_x={k:base0*(1-v/100) for k,v in P(red).items()},
                best_scaled_pct_reduction_same_layouts=P(redb),best_scaled_pct_reduction_all_layouts=P(redb_all)),red
s50,red50=stats(["50"]); sall,redall=stats(FV); out["how_much_k0_N200_b"]=dict(fv50=s50,pooled_all_mixes=sall)
fig,ax=plt.subplots(figsize=(12,7.5))
ax.hist(red50,bins=np.arange(0,12.01,0.5),color="#2a9d4a",edgecolor="k")
m=np.median(red50); b=red50.max()
ax.axvline(m,color="k",ls="--",lw=2); ax.text(m+0.15,ax.get_ylim()[1]*0.9,f"median: {m:.1f}% lower\n(237× → {base0*(1-m/100):.0f}×)",fontsize=13)
ax.axvline(b,color="#1d6e34",ls="-",lw=2); ax.text(b-0.15,ax.get_ylim()[1]*0.6,f"best: {b:.1f}% lower\n(237× → {base0*(1-b/100):.0f}×)",fontsize=13,ha="right")
ax.set_xlabel("how much lower the leftover got (percent below 237×)"); ax.set_ylabel("number of layouts")
ax.set_title(f"The {len(red50)} of 1,000 random layouts that helped (no clumps, 50% voids, 200 objects of 3×10⁴⁹–3×10⁵⁰ kg):\nhow much did each one help?",fontsize=15)
ax.grid(alpha=.3); fig.tight_layout(); fig.savefig("how_much_helped_k0.png",dpi=110); plt.close(fig)
out["ensemble_fitted"]=EN
json.dump(out,open("counts.json","w"),indent=1)
print(json.dumps(out["how_much_k0_N200_b"],indent=1))

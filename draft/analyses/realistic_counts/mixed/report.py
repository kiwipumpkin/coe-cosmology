import json, glob, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
plt.rcParams.update({"font.size":14,"axes.titlesize":15,"legend.fontsize":12})
from mix_core import *
from counts1x import analyse
from scan import CASES
# ---- collect scan results (keep the best = most converged per (case,X))
R={}
for f in ["log_scan_part1.txt","log_scan.txt","log_scan2.txt"]:
    for l in open(f):
        if not l.startswith("{"): continue
        r=json.loads(l); k=(r["case"],r["X"])
        if k not in R or r["left"]<R[k]["left"]: R[k]=r
json.dump([dict(v) for v in R.values()],open("scan_results.json","w"),indent=1)
# ---- 1x composition for threshold solutions (lambda=0 saved) + damped ones
C={}
for case,X,f in [("inside_mirror",0.065,"sol_inside_mirror_X0.065_F1.0_lam0.0_it1500.npz"),("both_nomirror",0.6,"sol_both_nomirror_X0.6_F1.0_lam0.0_it1500.npz")]:
    gs=CASES[case](); z=np.load(f); ds=[z[k] for k in sorted(z.files)]; r=analyse(gs,ds); r.update(X=X,lam=0.0,left=R[(case,X)]["left"]); C[(case,X)]=r
for case in ["inside_mirror","both_nomirror"]:
    if os.path.exists(f"counts1x_{case}.json"):
        for X,r in json.load(open(f"counts1x_{case}.json")).items(): C[(case,float(X))]=r
json.dump({f"{k[0]}|{k[1]}":v for k,v in C.items()},open("counts_at_1x.json","w"),indent=1)
for k in sorted(C): print(k,{a:(f"{b:.3g}" if isinstance(b,float) else b) for a,b in C[k].items()})
# ---- chart: best leftover vs allowed ratio deviation
COL={"inside_mirror":"#c0392b","both_nomirror":"#2471a3"}; LAB={"inside_mirror":"inside, mirror allowed","both_nomirror":"inside + outside, no mirror"}
fig,ax=plt.subplots(figsize=(10.5,7))
for case in COL:
    pts=sorted((X,max(v["left"],0.02)) for (c,X),v in R.items() if c==case)
    ax.plot([p[0]*100 for p in pts],[p[1] for p in pts],"o-",color=COL[case],lw=2.6,ms=7,label=LAB[case])
ax.axhline(1,color="green",lw=2); ax.text(1.05,1.12,"CMB flat target 1×",color="green",fontsize=12)
ax.axhline(237.3,color="0.5",ls="--",lw=1); ax.text(1.05,265,"start 237×",color="0.4",fontsize=11)
ax.axvspan(50,100,color="0.85",zorder=0); ax.text(70,140,"ratio 0.5×–2×\nof the mean",ha="center",fontsize=10,color="0.35")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.9,1000); ax.set_ylim(0.015,500)
ax.set_xlabel("allowed cluster:void ratio deviation from the cosmic mean (%)"); ax.set_ylabel("best leftover (× 10⁻⁵ floor)")
ax.set_title("Lobes built from a cluster + void mix: how far must the mix shift?")
ax.text(1.05,0.024,"points on the bottom line = fully cancelled (≈0.001×)",fontsize=10,color="0.4")
ax.legend(loc="center right",frameon=False); fig.tight_layout(); fig.savefig("leftover_vs_ratio_deviation.png",dpi=110); plt.close(fig)
# ---- map: inside threshold solution, per-lobe ratio
gs=CASES["inside_mirror"](); z=np.load("sol_inside_mirror_X0.065_F1.0_lam0.0_it1500.npz"); d=z[z.files[0]]; P=gs[0].P
a=d/DV_FRAC; b=np.zeros_like(d); q=1+a; lon=np.degrees(np.arctan2(P[:,1],P[:,0]))%360; r=np.linalg.norm(P,axis=1)
sel=np.abs(P[:,2])<0.26
fig,ax=plt.subplots(figsize=(10.5,10))
sc=ax.scatter(P[sel,0],P[sel,1],c=(q[sel]-1)*100,cmap="RdBu_r",vmin=-7,vmax=7,s=26,marker="s",linewidths=0)
ax.add_patch(Circle((0,0),14.3,fill=False,ls="--",color="k",lw=1.5)); ax.plot(0,0,"k+",ms=16,mew=2.5)
lobes=[]
for c,sg in [(15,1),(128,1),(195,1),(263,1),(75,-1),(162,-1),(228,-1),(315,-1)]:
    dl=np.abs((lon-c+180)%360-180); k=(dl<22)&(r>9)&(sg*d>0.003)
    if not k.any(): continue
    qm=float(np.average(q[k],weights=np.abs(d[k]))); vf=float(np.average(FV0*(1+b[k]),weights=np.abs(d[k])))
    nconv=float((a[k]*n_c0*gs[0].H**3).sum()); npres=float(((n_c0+n_v0)*gs[0].H**3)*k.sum())
    lobes.append(dict(lon=c,sign=sg,ratio=qm*NC_PER_CELL,ratio_rel=qm,void_fill=vf,converted=nconv,present=npres,mean_delta=float(d[k].mean())))
    aa=np.radians(c); txt=f"{qm*NC_PER_CELL:.2f} : 1\n({(qm-1)*100:+.1f}%)"
    ax.text(7.5*np.cos(aa),7.5*np.sin(aa),txt,ha="center",va="center",fontsize=12,color="#922b21" if sg>0 else "#1a5276",bbox=dict(facecolor="white",alpha=0.85,edgecolor="none",pad=1.5))
json.dump(lobes,open("lobes_mixed.json","w"),indent=1)
ax.set_xlim(-16,16); ax.set_ylim(-16,16); ax.set_aspect("equal"); ax.set_xlabel("Gpc"); ax.set_ylabel("Gpc")
cb=fig.colorbar(sc,ax=ax,shrink=0.7,pad=0.02); cb.set_label("cluster:void ratio vs cosmic mean 6.77:1 (%)")
ax.set_title("Best mixed solution (inside, mirror allowed, ratio within ±6.5%):\nclusters per void in each lobe (labels); slice through the ring plane",fontsize=13)
fig.tight_layout(); fig.savefig("map_mixed_ratio.png",dpi=110); plt.close(fig)
for l in lobes: print(l)

"""Charts + summary for the sequential (non-overlapping) scenarios. Plain linear charts."""
import os, json, glob, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
HERE=os.path.dirname(os.path.abspath(__file__)); os.chdir(HERE)
B0,B4=237.26,88.81
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False})

# ---------------- A
lf=json.load(open("local_fit.json")); best=[r for r in lf["runs"] if r["seed"]==lf["best_seed"]][0]
z=np.load("local_frozen.npz",allow_pickle=True); FROZEN=best["left"]
def pairs_traj(r):
    t=[p for p in r["traj"] if p["kind"] in ("start","cluster")]; return np.arange(len(t)),[p["left"] for p in t]
fig,ax=plt.subplots(figsize=(7,4.5))
for r in lf["runs"]:
    x,y=pairs_traj(r); ax.plot(x,y,color="#1f77b4",lw=2.5 if r is best else 1,alpha=1 if r is best else .35)
ax.axhline(B0,color="gray",ls="--"); ax.text(24,B0+4,"start 237×",ha="right",color="gray")
ax.axhline(B4,color="#d62728",ls="--"); ax.text(24,B4+4,"four clumps 89×",ha="right",color="#d62728")
ax.set_xlabel("void + cluster pairs placed near us (one at a time)"); ax.set_ylabel("leftover (× floor)")
ax.set_title(f"Scenario A step 1: near arrangement → {FROZEN:.0f}×"); ax.set_ylim(0,250); ax.set_xlim(0,24)
fig.tight_layout(); fig.savefig("a_local_fit.png",dpi=140); plt.close()

# layout picture (to scale, projected on the CMB-ring plane), shows nothing overlaps
def draw_layout(ax,X,r,kind,lim,title,slab=None):
    ax.add_patch(Circle((0,0),14.3,fill=False,ls="--",color="gray")); 
    for x,rr,k in zip(X,r,kind):
        if slab is not None and abs(x[2])>slab: continue
        ax.add_patch(Circle((x[0],x[1]),rr,color="#4a90d9" if k=="void" else "#d62728",alpha=.45 if k=="void" else .9,lw=0))
    ax.plot(0,0,"k+"); ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim); ax.set_aspect("equal"); ax.set_title(title,fontsize=10)
    ax.set_xlabel("Gpc"); 
fig,ax=plt.subplots(figsize=(6,6))
draw_layout(ax,z["X"],z["r"],z["kind"],32,"Frozen near arrangement (to scale, seen from above the CMB ring)\nblue = voids, red = clusters, dashed = 14.3 Gpc")
fig.tight_layout(); fig.savefig("a_local_layout.png",dpi=140); plt.close()

def gather(pat):
    rows=[]
    for f in sorted(glob.glob(pat)): rows+=json.load(open(f))["rows"]
    return rows
cfgs=[("empty (control)",None),("lattice\n(ordered)","outer_lattice_Rv1.5_s*.json"),("foam-like\n(grown)","outer_foam_Rv1.5_s*.json"),
      ("random\npacking","outer_random_Rv1.5_s1*.json"),("random,\nhalf-size voids","outer_random_Rv0.75_s*.json"),
      ("big pairs,\nscattered","outer_pairs_Rv0_s*.json")]
A={}
for lab,p in cfgs:
    if p is None: A[lab]=dict(n=1,left=np.array([FROZEN]),outer=np.array([0.0]),nolocal=np.array([B0]),units=0,vf=0); continue
    rows=gather(p)
    if not rows: continue
    A[lab]=dict(n=len(rows),left=np.array([q["left"] for q in rows]),outer=np.array([q["outer_rms"] for q in rows]),
                nolocal=np.array([q["left_without_local"] for q in rows]),units=np.median([q["units"] for q in rows]),vf=np.median([q["void_vol_frac"] for q in rows]),
                gap=min(q["min_gap"] for q in rows))
labs=list(A)
fig,ax=plt.subplots(figsize=(9,4.8))
med=[np.median(A[l]["left"]) for l in labs]; lo=[np.percentile(A[l]["left"],10) for l in labs]; hi=[np.percentile(A[l]["left"],90) for l in labs]
cols=["#9e9e9e","#ff7f0e","#2ca02c","#1f77b4","#17becf","#9467bd"][:len(labs)]
ax.bar(range(len(labs)),med,color=cols,width=.6)
ax.errorbar(range(len(labs)),med,yerr=[np.array(med)-lo,np.array(hi)-med],fmt="none",ecolor="k",capsize=4,lw=1)
for i,m in enumerate(med): ax.text(i,hi[i]+3,f"{m:.1f}×",ha="center",fontsize=10)
ax.axhline(FROZEN,color="k",ls=":"); ax.axhline(B4,color="#d62728",ls="--"); ax.text(len(labs)-.5,B4+3,"four clumps 89×",ha="right",color="#d62728",fontsize=9)
ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs,fontsize=9); ax.set_ylabel("leftover (× floor)")
ax.set_title(f"Scenario A step 2: near arrangement frozen at {FROZEN:.1f}×, then the outer shell filled\n(bar = median, line = middle 80% of layouts)",fontsize=11)
ax.set_ylim(0,max(hi)*1.12); fig.tight_layout(); fig.savefig("a_outer_effect.png",dpi=140); plt.close()

fig,ax=plt.subplots(figsize=(8,4.2)); L2=[l for l in labs if A[l]["n"]>1]
h=[int(np.sum(A[l]["left"]<FROZEN-0.01)) for l in L2]; u=[int(np.sum(A[l]["left"]>FROZEN+0.01)) for l in L2]; s=[A[l]["n"]-a-b for l,a,b in zip(L2,h,u)]
x=np.arange(len(L2)); ax.bar(x-.27,h,.27,color="#2ca02c",label="helped"); ax.bar(x,u,.27,color="#d62728",label="hurt"); ax.bar(x+.27,s,.27,color="#bdbdbd",label="no change (<0.01×)")
for i in range(len(L2)):
    for dx,v in ((-.27,h[i]),(0,u[i]),(.27,s[i])): ax.text(i+dx,v+1,str(v),ha="center",fontsize=9)
ax.set_xticks(x); ax.set_xticklabels(L2,fontsize=9); ax.set_ylabel("layouts"); ax.legend(frameon=False)
ax.set_title(f"Outer-shell fills: helped vs hurt compared with the frozen {FROZEN:.1f}×"); fig.tight_layout(); fig.savefig("a_helped_hurt.png",dpi=140); plt.close()

# ---------------- B
B={}
for f in sorted(glob.glob("b_*_Rv*_s*.json")):
    d=json.load(open(f)); B.setdefault((d["rule"],d["Rv"]),[]).append(d)
RL={"foam_rim":"foam: voids grow next to voids, clusters at junctions","alternate":"alternate: void, then one heavy cluster across a neighbour",
    "coe_bias":"CoE-biased: voids away from CoE, clusters facing it","random":"random: same objects dropped anywhere (control)",
    "cmb_aware":"CMB-aware: each spot chosen to lower the leftover"}
CB={"foam_rim":"#2ca02c","alternate":"#ff7f0e","coe_bias":"#9467bd","random":"#7f7f7f","cmb_aware":"#1f77b4"}
fig,ax=plt.subplots(figsize=(9,5))
for rule in RL:
    for i,d in enumerate(B.get((rule,1.5),[])):
        ax.plot(d["traj_void"],d["traj_left"],color=CB[rule],lw=1.6 if i==0 else .8,alpha=1 if i==0 else .4,label=RL[rule] if i==0 else None)
ax.axhline(B0,color="gray",ls="--"); ax.axhline(B4,color="#d62728",ls="--"); ax.text(50,B4-40,"four clumps 89×",color="#d62728",fontsize=9)
ax.set_xlabel("voids placed so far (each followed by its cluster mass)"); ax.set_ylabel("leftover (× floor)")
ax.set_title("Scenario B: leftover as the shell is grown one object at a time\n(void radius 1.5 Gpc; dashed grey = start 237×)",fontsize=11)
ax.set_ylim(0,None); ax.legend(frameon=False,fontsize=8,loc="upper left"); fig.tight_layout(); fig.savefig("b_growth.png",dpi=140); plt.close()

# zoom: first 300 voids
fig,ax=plt.subplots(figsize=(9,5))
for rule in RL:
    for i,d in enumerate(B.get((rule,1.5),[])):
        v=np.array(d["traj_void"]); k=v<=300
        ax.plot(v[k],np.array(d["traj_left"])[k],color=CB[rule],lw=1.6 if i==0 else .8,alpha=1 if i==0 else .4,label=RL[rule] if i==0 else None)
ax.axhline(B0,color="gray",ls="--"); ax.axhline(B4,color="#d62728",ls="--")
ax.set_xlabel("voids placed so far"); ax.set_ylabel("leftover (× floor)"); ax.set_title("Scenario B, first 300 voids (zoom)",fontsize=11)
ax.set_ylim(0,None); ax.legend(frameon=False,fontsize=8); fig.tight_layout(); fig.savefig("b_growth_first300.png",dpi=140); plt.close()

fig,ax=plt.subplots(figsize=(9,4.5)); rules=[r for r in RL if (r,1.5) in B]
fin=[np.median([d["final_left"] for d in B[(r,1.5)]]) for r in rules]; mn=[np.median([d["min_left"] for d in B[(r,1.5)]]) for r in rules]
x=np.arange(len(rules)); ax.bar(x-.2,fin,.4,color=[CB[r] for r in rules],label="when the shell is full")
ax.bar(x+.2,mn,.4,color=[CB[r] for r in rules],alpha=.45,hatch="//",label="lowest point during growth")
for i in range(len(rules)): ax.text(i-.2,fin[i]+8,f"{fin[i]:.0f}×",ha="center",fontsize=9); ax.text(i+.2,mn[i]+8,f"{mn[i]:.0f}×",ha="center",fontsize=9)
ax.axhline(B0,color="gray",ls="--"); ax.axhline(B4,color="#d62728",ls="--")
ax.set_xticks(x); ax.set_xticklabels([RL[r].split(":")[0] for r in rules]); ax.set_ylabel("leftover (× floor), median of seeds")
ax.legend(frameon=False,fontsize=9); ax.set_title("Scenario B: where each growth rule ends up"); fig.tight_layout(); fig.savefig("b_final_bars.png",dpi=140); plt.close()

# example B layouts: thin slab through the ring plane
from scenario_b import Grow
fig,axs=plt.subplots(1,2,figsize=(11,5.6))
for ax,rule in zip(axs,("foam_rim","random")):
    g=Grow(rule,7,1.5); g.max_voids=2500; g.run()
    Xa,ra,Ma,ka=g.pk.arrays()
    draw_layout(ax,Xa,ra,ka,44,f"{rule}: first 2,500 voids, slice |z|<1.5 Gpc (to scale)",slab=1.5)
fig.tight_layout(); fig.savefig("b_example_slices.png",dpi=120); plt.close()

json.dump(dict(frozen=FROZEN,local_runs=[dict(seed=r["seed"],left=r["left"],min_gap=r["min_gap"],n_obj=r["n_obj"],final_amps=r["final_amps"]) for r in lf["runs"]],
               start_amps=lf["start_amps"],four_clump_amps=lf["four_clump_amps"],
               outer={l.replace("\n"," "):dict(n=A[l]["n"],median=float(np.median(A[l]["left"])),p10=float(np.percentile(A[l]["left"],10)),p90=float(np.percentile(A[l]["left"],90)),
                     min=float(A[l]["left"].min()),max=float(A[l]["left"].max()),outer_rms_median=float(np.median(A[l]["outer"])),
                     same_fill_without_local_median=float(np.median(A[l]["nolocal"])),helped=int(np.sum(A[l]["left"]<FROZEN-0.01)),hurt=int(np.sum(A[l]["left"]>FROZEN+0.01)),
                     units=float(A[l]["units"]),void_vol_frac=float(A[l]["vf"]),min_gap=float(A[l].get("gap",np.inf))) for l in labs},
               B={f"{k[0]}_Rv{k[1]:g}":dict(n=len(v),final=[d["final_left"] for d in v],min=[d["min_left"] for d in v],web_rms=[d["web_rms"] for d in v],
                   n_void=[d["n_void"] for d in v],n_clu=[d["n_clu"] for d in v],void_vol_frac=[d["void_vol_frac"] for d in v],min_gap=[d["min_gap_gpc"] for d in v],
                   net_kg=[d["net_kg"] for d in v]) for k,v in B.items()}),open("results.json","w"),indent=1)
print("ok")

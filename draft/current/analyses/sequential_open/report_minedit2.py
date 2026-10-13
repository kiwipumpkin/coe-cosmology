"""Report for minedit2 (move / rescale / delete existing bodies)."""
import os,json,pickle,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from common import *
plt.rcParams.update({"font.size":7.5})
RUNS=[("b_random_s1_tidalcore","random s1","#95a5a6"),("b_random_s2_tidalcore","random s2","#717d7e"),("b_random_s3_tidalcore","random s3","#2c3e50"),
      ("b_foam_rim_s1_tidalcore","foam s1","#85c1e9"),("b_foam_rim_s2_tidalcore","foam s2","#3498db"),("b_foam_rim_s3_tidalcore","foam s3","#1b4f72")]
import glob as _g,re as _re
_base=list(RUNS); RUNS=list(_base)
for t,nm,c in _base:
    Ds=sorted({float(_re.search(r"_clr([0-9.]+)",p).group(1)) for p in _g.glob(f"{OUT}/minedit2/{t}_clr*.json")+_g.glob(f"{OUT}/ckpt/minedit2_{t}_clr*.pkl") if "pruned" not in p})
    RUNS+=[(f"{t}_clr{D:g}",f"{nm}, ring zone cleared to {D:g} Gpc first",c) for D in Ds]
BASE=lambda t:t.split("_clr")[0]
from check_layout import check
def dring(X): X=np.atleast_2d(X); return np.sqrt((np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2)
D={}
for t,nm,c in RUNS:
    fj=f"{OUT}/minedit2/{t}.json"; ck=f"{OUT}/ckpt/minedit2_{t}.pkl"
    if os.path.exists(fj): D[t]=json.load(open(fj)); D[t]["final"]=True
    elif os.path.exists(ck): D[t]=dict(curve=pickle.load(open(ck,"rb"))["curve"],final=False)
runs=[r for r in RUNS if r[0] in D]; fin=[r for r in runs if D[r[0]]["final"]]
V1={}
for t,_,_ in RUNS:
    p=f"{OUT}/minedit/{BASE(t)}_nr.json"
    if os.path.exists(p): V1[t]=json.load(open(p))["curve"]
fig,(ax,ax2)=plt.subplots(2,1,figsize=(7,9.5),gridspec_kw=dict(height_ratios=[2,1]))
for t,nm,c in _base:
    if t in D:
        cv=np.array(D[t]["curve"]); ax.plot(np.maximum(cv[:,0],0.8),cv[:,1],"-",color=c,lw=1.5,label=nm+(": edits only (move/resize/delete)" if D[t]["final"] else " (running)"))
    pts=sorted([(D[k]["n_edited"],D[k]["left"],D[k]["clear_D"]) for k,_,_ in runs if k.startswith(t+"_clr") and D[k]["final"]])
    if pts:
        p=np.array(pts); ax.plot(p[:,0],p[:,1],"--o",color=c,lw=1.0,ms=4,label=nm+": clear ring zone to D, then edit")
        ax2.plot(p[:,0],p[:,1],"--o",color=c,lw=1.0,ms=4)
        for x_,y_,d_ in pts: ax2.annotate(f"{d_:g}",(x_,y_),textcoords="offset points",xytext=(3,3),fontsize=6,color=c)
ax.axhline(1,color="green",lw=1.5); ax.set_xscale("log"); ax.set_yscale("log")
ax2.axhline(1,color="green",lw=1.5); ax2.set_xscale("log"); ax2.set_ylim(0,2.5); ax2.set_xlabel("existing bodies edited, incl. the cleared ones (log); labels = cleared zone D (Gpc)"); ax2.set_ylabel("leftover (× 1×, linear)"); ax2.set_title("zoom: clear-then-edit runs (one point per run)",fontsize=8)
ax.set_xlabel("existing bodies edited (moved, resized or deleted; log, 0 at 0.8); numbers = cleared zone D (Gpc)"); ax.set_ylabel("leftover (× 1×)")
ax.set_title("Minimal edit with deletion allowed (joint refit of the selected set)",fontsize=8.5)
hd=[l for l in ax.get_lines() if not l.get_label().startswith("_")]+[Line2D([],[],color="0.4",ls=":",lw=0.9,label="earlier move/resize-only pass (had an overlap bug)"),Line2D([],[],color="green",lw=1.5,label="1× (CMB level)")]
fig.legend(handles=hd,loc="upper center",ncol=2,frameon=False,fontsize=6.5,bbox_to_anchor=(0.5,0.0)); fig.tight_layout(); fig.savefig(f"{OUT}/minedit2_curve.png",dpi=200,bbox_inches="tight"); plt.close(fig)
for grp,fn_ in (([r for r in fin if '_clr' not in r[0]],'minedit2_changes.png'),([r for r in fin if '_clr' in r[0]],'minedit2_changes_cleared.png')):
  if not grp: continue
  fin_=grp
  if True:
    fig,axs=plt.subplots(3,2,figsize=(7.4,9.4))
    w=0.8/len(fin_)
    for k,(t,nm,c) in enumerate(fin_):
        d=D[t]; de=np.array(d["deleted"]); mo=np.array(d["moved"]); rs=np.array(d["resized"]); nr=np.array(d["near_ring"]); iv=np.array(d["is_void"])
        cats=[(de&iv).sum(),(de&~iv).sum(),(mo&iv).sum(),(mo&~iv).sum(),(rs&iv).sum(),(rs&~iv).sum()]
        axs[0,0].bar(np.arange(6)+k*w-0.4+w/2,cats,width=w,color=c)
        axs[0,1].bar(np.arange(3)+k*w-0.4+w/2,[(de&nr).sum(),(mo&nr).sum(),(rs&nr&~mo).sum()],width=w,color=c)
        mv=np.array(d["move"])[mo]; axs[1,0].hist(mv,bins=np.linspace(0,15,31),histtype="step",color=c,lw=1.2)
        a=np.array(d["a_new"])[rs]; axs[1,1].hist(np.log10(a),bins=np.linspace(-1.05,1.05,22),histtype="step",color=c,lw=1.2)
        keep=~de; axs[2,0].hist(np.array(d["dist_us_old"])[de],bins=np.linspace(0,30,31),histtype="step",color=c,lw=1.2,ls="--")
        axs[2,0].hist(np.array(d["dist_us_new"])[keep],bins=np.linspace(0,30,31),histtype="step",color=c,lw=1.2)
        axs[2,1].hist(np.array(d["dist_ring_old"])[de],bins=np.linspace(0,25,26),histtype="step",color=c,lw=1.2,ls="--")
        axs[2,1].hist(np.array(d["dist_ring_new"])[keep],bins=np.linspace(0,25,26),histtype="step",color=c,lw=1.2)
    axs[0,0].set_xticks(range(6)); axs[0,0].set_xticklabels(["deleted\nvoids","deleted\nclusters","moved\nvoids","moved\nclusters","resized\nvoids","resized\nclusters"],fontsize=6.5)
    axs[0,1].set_xticks(range(3)); axs[0,1].set_xticklabels(["near-ring\ndeleted","near-ring\nmoved out","near-ring\nresized only"],fontsize=6.5)
    axs[1,0].set_xlabel("move distance of moved bodies (Gpc)"); axs[1,1].set_xlabel("log10 mass scale of resized bodies (kept ones)")
    axs[2,0].set_xlabel("distance from us (Gpc): kept edits (solid, new), deleted (dashed, old)"); axs[2,1].set_xlabel("distance from the ring (Gpc): kept (solid), deleted (dashed)")
    f0=np.load(f"{OUT}/final_{BASE(fin_[0][0])}.npz"); Xa=f0["Xp"]; sc=np.mean([D[t]["n_edited"] for t,_,_ in fin_])/len(Xa)
    axs[2,0].hist(np.linalg.norm(Xa,axis=1),bins=np.linspace(0,30,31),color="green",alpha=0.2,weights=np.full(len(Xa),sc))
    axs[2,1].hist(dring(Xa),bins=np.linspace(0,25,26),color="green",alpha=0.2,weights=np.full(len(Xa),sc)); axs[2,1].axvline(3,color="red",ls=":",lw=0.8)
    for a_ in axs.ravel(): a_.set_ylabel("bodies")
    hd=[plt.Rectangle((0,0),1,1,fc=c,label=nm) for t,nm,c in fin_]+[plt.Rectangle((0,0),1,1,fc="green",alpha=0.2,label="all web bodies (scaled to the edit count)")]
    fig.legend(handles=hd,loc="upper center",ncol=4,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.0)); fig.tight_layout(); fig.savefig(f"{OUT}/{fn_}",dpi=200,bbox_inches="tight"); plt.close(fig)
    if fn_!='minedit2_changes.png': continue
    if not fin: continue
    from edit_before_after import figure as _ba
    eo=sorted([r for r in fin if "_clr" not in r[0]],key=lambda r:D[r[0]]["left"])[:1]
    d2=sorted([r for r in fin if r[0].endswith("_clr2")],key=lambda r:D[r[0]]["left"])[:1]
    rc=sorted([r for r in fin if "_clr" in r[0] and D[r[0]]["reached"] and "random" in r[0]],key=lambda r:D[r[0]]["n_edited"])[:1]
    fc=sorted([r for r in fin if "_clr" in r[0] and D[r[0]]["reached"] and "foam" in r[0]],key=lambda r:D[r[0]]["n_edited"])[:1]
    _ba([(f"{OUT}/minedit2/{t}.json",nm,c) for t,nm,c in eo+d2+rc+fc],f"{OUT}/minedit2_side_views.png")
S={}
for t,_,_ in runs:
    d=D[t]
    if not d["final"]: S[t]=dict(running=True,last=d["curve"][-1]); continue
    a=d["after"]; ch=np.array(d["moved"])
    S[t]=dict(n_bodies=d["n_bodies"],edited=d["n_edited"],frac=round(d["frac_edited"],4),deleted=d["n_deleted"],moved=d["n_moved"],resized=d["n_resized"],near_ring=d["n_near_ring"],near_ring_edited=d["n_near_ring_edited"],
      before=round(d["before"]["total"],1),left=round(d["left"],2),reached=d["reached"],m1_12=round(float(np.hypot(a["rms_m1_5"],a["rms_m6_12"])),2),m13=round(a["rms_13up"],2),sym=round(a["sym"],2),
      legality=check(f"{OUT}/minedit2/{t}.json",BASE(t)),cleared=d.get("n_cleared",0),
      median_move=round(float(np.median(np.array(d["move"])[ch])),2) if ch.any() else 0,delete_all_near_ring_only=round(d["delete_all_near_ring"]["left"],1),net_mass_change=d["net_mass_change"])
json.dump(S,open(f"{OUT}/minedit2/summary.json","w"),indent=1)
for k,v in S.items(): print(k,v)

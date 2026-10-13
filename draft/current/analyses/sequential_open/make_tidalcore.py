"""Four growth runs, 30 Gpc fill INCLUDING the part of the CoE core inside 30 Gpc (rest of core empty). Curves vs old and vs core-empty tidal runs + 4 final universes."""
import os,json,glob,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from common import *
plt.rcParams.update({"font.size":8,"savefig.dpi":300,"pdf.fonttype":42})
RUNS=[("seqB_random","b","random",2,"#7f8c8d","B: growth, random"),("seqB_foam_rim","b","foam_rim",3,"#5dade2","B: growth, foam rule"),
      ("seqB_cmb_aware","b","cmb_aware",1,"#7d3c98","B: growth, CMB-aware"),("seqB_shallow","shallow","cmb_aware",2,"#b48fc9","B: growth, CMB-aware, shallow voids")]
VAR=[("_tidal","30 Gpc fill, core empty",(0,(3,1.5))),("_tidalcore","30 Gpc fill incl. core part inside 30 Gpc",'-')]
def old_curve(kind,rule,seed):
    if kind=="b": d=json.load(open(f"{SEQ}/b_{rule}_Rv1.5_s{seed}.json")); return d["traj_obj"],d["traj_left"],d["final_left"],d["min_left"],d["n_obj"]
    d=json.load(open(f"{SEQ}/shallow/grow_{rule}_s{seed}.json")); return [5*i+1 for i in range(len(d["traj"]))],d["traj"],d["final"],d["min"],d["n_void"]
def new(kind,rule,seed,suf):
    fj=f"{OUT}/{'b' if kind=='b' else 'shallow'}_{rule}_s{seed}{suf}.json"
    if not os.path.exists(fj): return None
    d=json.load(open(fj)); f=np.load(f"{OUT}/final_{'b' if kind=='b' else 'shallow'}_{rule}_s{seed}{suf}.npz")
    if kind=="b":
        v=f["kp"]=="void"; L=dict(Xv=f["Xp"][v],rv=f["rp"][v],Xc=f["Xp"][~v],n=f"{v.sum():,} voids + {(~v).sum():,} clusters",wa=d["where_all"],X=f["X"],M=f["M"])
        L.update(final=d["final_left"],best=d["min_left"],at=d["argmin"],nobj=d["n_obj"],vf=d["void_vol_frac"],gap=d["min_gap_gpc"])
    else:
        L=dict(Xv=f["X"],rv=np.full(len(f["X"]),1.5),Xc=np.zeros((0,3)),n=f"{len(f['X']):,} shallow voids",wa=d["where_voids"],X=f["X"],M=None)
        L.update(final=d["final"],best=d["min"],at=d["argmin"],nobj=d["n_void"],vf=d["void_vol_frac"],gap=d["min_gap"],cap=d["cap"],hit_cap=d["hit_cap"])
    L.update(w=f["w"],x=d["traj_obj"],y=d["traj_left"]); return L
D={}; OLD={}
for key,kind,rule,seed,col,nm in RUNS:
    OLD[key]=old_curve(kind,rule,seed)
    for suf,_,_ in VAR: D[(key,suf)]=new(kind,rule,seed,suf)
# ---- groups: what each region's bodies contribute (final layouts)
def deep_k(X,M):
    W=np.zeros(360)
    for i in range(0,len(X),3000): W+=M[i:i+3000]@seq_core.unit_pattern_fast(X[i:i+3000])
    return W
import shallow_kern as sk
G={}
for (key,suf),L in D.items():
    if L is None: continue
    X=L["X"]; s=np.linalg.norm(X,axis=1); dc=np.linalg.norm(X-C_G,axis=1)
    grp={"inside_us":s<RMIN,"core":(dc<A_IN_G)&(s>=RMIN),"shell":(dc>=A_IN_G)&(s>=RMIN)}
    K=(lambda m: deep_k(X[m],L["M"][m])) if L["M"] is not None else (lambda m: sk.pat(X[m]))
    W={g:(K(m) if m.any() else np.zeros(360)) for g,m in grp.items()}; tot=sum(W.values())
    G[f"{key}{suf}"]={g:dict(n=int(grp[g].sum()),own=float(rms(W[g])/FLOOR),without=float(rms(R-tot+W[g])/FLOOR)) for g in grp}
    G[f"{key}{suf}"]["total"]=float(rms(R-tot)/FLOOR)
json.dump(G,open(f"{OUT}/groups_tidalcore.json","w"),indent=1)
# ---- curves: one panel per run
fig,axs=plt.subplots(2,2,figsize=(7.1,8.4),sharex=True,sharey=True)
for ax,(key,kind,rule,seed,col,nm) in zip(axs.ravel(),RUNS):
    x,y=OLD[key][:2]; ax.plot(x,y,color="0.45",lw=1.0,ls=(0,(4,2)),label="old (14.3–42.9 Gpc, core blocked)")
    for suf,lab,ls in VAR:
        L=D[(key,suf)]
        if L is not None: ax.plot(L["x"],L["y"],color=col,lw=1.5,ls=ls,label=lab)
    Lb=new(kind,rule,seed,"")
    if Lb is not None: ax.plot(Lb["x"],Lb["y"],color="0.6",lw=0.9,ls=(0,(1,1.5)))
    ax.axhline(1,color="green",lw=1.5); ax.axhline(237.3,color="0.5",ls=":",lw=0.8); ax.axhline(88.8,color="0.6",ls=":",lw=0.8)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1,2e5); ax.set_ylim(0.5,5000); ax.set_title(nm,color=col,fontweight="bold",fontsize=8)
    ax.text(1.2,1.1,"1×",color="green",fontsize=6.5,va="bottom"); ax.text(1.2,250,"237×",color="0.4",fontsize=6.5,va="bottom"); ax.text(1.2,93,"89×",color="0.45",fontsize=6.5,va="bottom")
for ax in axs[1]: ax.set_xlabel("bodies placed")
for ax in axs[:,0]: ax.set_ylabel(r"leftover ($\times$ floor)")
hd=[Line2D([],[],color="0.45",lw=1.0,ls=(0,(4,2)),label="old run (14.3–42.9 Gpc, core blocked)")]+[Line2D([],[],color="0.6",lw=0.9,ls=(0,(1,1.5)),label="earlier: everything filled incl. core (superseded)")]+[Line2D([],[],color="0.2",lw=1.3,ls=ls,label=lab) for suf,lab,ls in VAR]+[Line2D([],[],color="green",lw=1.5,label="CMB flat, 1×")]
fig.legend(handles=hd,loc="lower center",ncol=2,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.0)); fig.tight_layout(rect=(0,0.07,1,1))
fig.savefig(f"{OUT}/fig7_growth_tidalcore.png",dpi=200,bbox_inches="tight"); fig.savefig(f"{OUT}/fig7_growth_tidalcore.pdf",bbox_inches="tight"); plt.close(fig)
# ---- universes: rows = runs (largest 'both' leftover first), columns = variants
rows=sorted(RUNS,key=lambda r:-(D[(r[0],"_tidalcore")]["final"] if D[(r[0],"_tidalcore")] else 0))
fig=plt.figure(figsize=(10.5,11)); GO=fig.add_gridspec(2,2,hspace=0.32,wspace=0.18,left=0.06,right=0.98,top=0.96,bottom=0.08)
for i,(key,kind,rule,seed,col,nm) in enumerate(rows):
    for j,(suf,lab,_) in enumerate(VAR[1:]):
        L=D[(key,suf)]; sub=GO[i//2,i%2].subgridspec(2,1,height_ratios=[1.6,0.9],hspace=0.55); ax=fig.add_subplot(sub[0])
        if L is None: ax.set_title(f"{nm}\n{lab}: still running",fontsize=7); ax.axis("off"); fig.add_subplot(sub[1]).axis("off"); continue
        core_open=False; ax.add_patch(Circle((0,0),30.0,fc="#fadbd8",ec="none",zorder=0,alpha=0.0))
        cp=Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="none",zorder=1); ax.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=ax.transData))
        Xv,rv,Xc=L["Xv"],L["rv"],L["Xc"]; sv=(np.abs(Xv[:,2])<1.5)&(np.linalg.norm(Xv,axis=1)<=30); sc=(np.abs(Xc[:,2])<1.5)&(np.linalg.norm(Xc,axis=1)<=30)
        ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(Xv[sv],rv[sv])],fc="#aed6f1",ec="#2471a3",lw=0.2,zorder=2))
        if len(Xc): ax.scatter(Xc[sc,0],Xc[sc,1],s=0.8,c="#c0392b",lw=0,zorder=3)
        ax.add_patch(Circle((0,0),RMIN,fill=suf!="_core" and False,ec="#1f4e79",lw=1.3,zorder=4)); ax.plot(0,0,"k+",ms=5,mew=1,zorder=5)
        frame30(ax); ax.tick_params(labelsize=5.5)
        for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(2.0)
        ax.set_title(f"{nm} — {lab}\n{L['n']} (slab |z|<1.5 Gpc)",color=col,fontsize=6.5,fontweight="bold")
        r_=fig.add_subplot(sub[1]); res=(R-L["w"])/FLOOR
        r_.axhspan(-1,1,color="green",alpha=0.25,lw=0); r_.axhline(1,color="green",lw=0.7); r_.axhline(-1,color="green",lw=0.7); r_.axhline(0,color="0.6",lw=0.5)
        r_.plot(np.arange(360),res,color=col,lw=1.0); lim=1.3*max(np.abs(res).max(),2); r_.set_yscale('symlog',linthresh=1,linscale=0.6); r_.set_ylim(-lim,lim); r_.set_xlim(0,360); r_.set_xticks([0,90,180,270,360]); r_.tick_params(labelsize=5.5)
        r_.set_title(f"leftover {L['final']:,.1f}× (lowest {L['best']:,.1f}×)",fontsize=6.5,color=col,fontweight="bold"); r_.set_xlabel("ℓ (°) around the ring",fontsize=5.8,labelpad=1); r_.set_ylabel("× $10^{-5}$\n(linear in ±1, log beyond)",fontsize=5.5)
hd=[Line2D([],[],marker="o",ls="",color="#c0392b",ms=3.5,label="cluster"),Line2D([],[],marker="o",ls="",mfc="#aed6f1",mec="#2471a3",ms=7,label="void (to scale)"),
    Line2D([],[],color="#1f4e79",lw=1.4,label="CMB ring / our sphere (14.3 Gpc)"),
    plt.Rectangle((0,0),1,1,fc="#fadbd8",label="part of the core within 30 Gpc: filled"),Line2D([],[],color="0.4",lw=0.8,ls=":",label="fill limit: 30 Gpc from us (tidal depth)"),Line2D([],[],color="#a93226",marker=">",lw=1.2,label="direction to the CoE (core beyond 30 Gpc not shown)"),
    plt.Rectangle((0,0),1,1,fc="green",alpha=0.25,ec="none",label="±1× (CMB level)")]
fig.legend(handles=hd,loc="lower center",ncol=4,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.0))
fig.savefig(f"{OUT}/fig8_growth_tidalcore.png",dpi=170,bbox_inches="tight"); fig.savefig(f"{OUT}/fig8_growth_tidalcore.pdf",bbox_inches="tight"); plt.close(fig)
T={}
for key,kind,rule,seed,col,nm in RUNS:
    T[key]=dict(old=dict(final=OLD[key][2],best=OLD[key][3],n=OLD[key][4]))
    for suf,lab,_ in VAR:
        L=D[(key,suf)]
        if L: T[key][suf]=dict(final=L["final"],best=L["best"],at=L["at"],n=L["nobj"],vf=L["vf"],where=L["wa"],groups=G.get(f"{key}{suf}"))
json.dump(T,open(f"{OUT}/results_tidalcore.json","w"),indent=1)
for k,v in T.items():
    print(k,"old",{q:round(x,1) for q,x in v["old"].items()})
    for vv in ("_tidal","_tidalcore"):
        if vv in v: x=v[vv]; print("   ",vv,"final %.1f best %.1f at %d n %d vf %.3f"%(x["final"],x["best"],x["at"],x["n"],x["vf"]),x["where"],{g:(y["n"],round(y["own"],1),round(y["without"],1)) for g,y in x["groups"].items() if g!="total"})

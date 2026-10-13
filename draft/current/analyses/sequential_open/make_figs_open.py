"""Fig-7-style curves (open vs old) and Fig-8-style six final universes for the open sequential runs."""
import os,sys,json,glob,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from common import *
plt.rcParams.update({"font.size":8,"savefig.dpi":300,"pdf.fonttype":42})
S_MAX=42.9
STY={"seqA_deep":"#c0392b","seqA_shallow":"#e59866","seqB_cmb_aware":"#7d3c98","seqB_foam_rim":"#5dade2","seqB_random":"#7f8c8d","seqB_shallow":"#d7bde2"}
NAME={"seqA_deep":"A: near fix, deep voids + clusters","seqA_shallow":"A: near fix, shallow voids","seqB_cmb_aware":"B: growth, CMB-aware",
      "seqB_foam_rim":"B: growth, foam rule","seqB_random":"B: growth, random","seqB_shallow":"B: growth, CMB-aware, shallow voids"}
SEEDS={"seqB_cmb_aware":("cmb_aware",1),"seqB_foam_rim":("foam_rim",3),"seqB_random":("random",2)}
# ---------- old curves (exactly as the brief's collect.py)
OLD={}
lf=json.load(open(f"{SEQ}/local_fit.json")); bs=min(lf["runs"],key=lambda r:r["left"]); OLD["seqA_deep"]={i:t["left"] for i,t in enumerate(bs["traj"]) if i>0}
nf=json.load(open(f"{SEQ}/shallow/near_fit.json"))["best"]; OLD["seqA_shallow"]={i:v for i,v in enumerate(nf["traj"]) if i>0}
for k,(rule,_) in SEEDS.items():
    runs=[json.load(open(f)) for f in glob.glob(f"{SEQ}/b_{rule}_Rv1.5_s*.json")]; d=min(runs,key=lambda r:r["min_left"]); OLD[k]=dict(zip(d["traj_obj"],d["traj_left"]))
g=[json.load(open(f)) for f in glob.glob(f"{SEQ}/shallow/grow_cmb_aware_s*.json")]; d=min(g,key=lambda r:r["min"]); OLD["seqB_shallow"]={5*i:v for i,v in enumerate(d["traj"]) if i>0}
# ---------- new
NEW={}; E={}; INFO={}
a=json.load(open(f"{OUT}/deepA_near.json")); b=min(a["runs"],key=lambda r:r["left"]); NEW["seqA_deep"]={i:v for i,v in enumerate(b["traj"]) if i>0}
z=np.load(f"{OUT}/deepA_best.npz",allow_pickle=True); k=z["kind"].astype(str)
E["seqA_deep"]=dict(Xv=z["X"][k=="void"],rv=z["r"][k=="void"],Xc=z["X"][k=="clu"],w=z["w"],n=f"{(k=='void').sum()} deep voids + {(k=='clu').sum()} clusters",slab=False)
INFO["seqA_deep"]=dict(final=b["left"],best=min(b["traj"][1:]),n=b["n_obj"],where=b["where"],seeds=[round(r["left"],2) for r in a["runs"]])
s=json.load(open(f"{OUT}/shallowA_near.json")); bb=s["best"]; NEW["seqA_shallow"]={i:v for i,v in enumerate(bb["traj"]) if i>0}
E["seqA_shallow"]=dict(Xv=np.array(bb["X"]),rv=np.array(bb["Rv"]),Xc=np.zeros((0,3)),w=np.array(bb["w"]),n=f"{bb['n']} shallow voids",slab=False)
INFO["seqA_shallow"]=dict(final=bb["left"],best=min(bb["traj"][1:]),n=bb["n"],where=where(np.array(bb["X"])),seeds=[round(r["left"],2) for r in s["runs"]])
for kk,(rule,seed) in SEEDS.items():
    d=json.load(open(f"{OUT}/b_{rule}_s{seed}.json")); NEW[kk]=dict(zip(d["traj_obj"],d["traj_left"]))
    f=np.load(f"{OUT}/final_b_{rule}_s{seed}.npz"); v=f["kp"]=="void"
    E[kk]=dict(Xv=f["Xp"][v],rv=f["rp"][v],Xc=f["Xp"][~v],w=f["w"],n=f"{v.sum():,} voids + {(~v).sum():,} clusters",slab=True)
    INFO[kk]=dict(final=d["final_left"],best=d["min_left"],at=d["argmin"],n=d["n_obj"],where=d["where_all"],wv=d["where_voids"],wc=d["where_clusters"],vf=d["void_vol_frac"],gap=d["min_gap_gpc"])
d=json.load(open(f"{OUT}/shallow_cmb_aware_s2.json")); NEW["seqB_shallow"]=dict(zip(d["traj_obj"],d["traj_left"]))
f=np.load(f"{OUT}/final_shallow_cmb_aware_s2.npz"); E["seqB_shallow"]=dict(Xv=f["X"],rv=np.full(len(f["X"]),1.5),Xc=np.zeros((0,3)),w=f["w"],n=f"{len(f['X']):,} shallow voids",slab=True)
INFO["seqB_shallow"]=dict(final=d["final"],best=d["min"],at=d["argmin"],n=d["n_void"],where=d["where_voids"],vf=d["void_vol_frac"],cap=d["cap"],hit_cap=d["hit_cap"],gap=d["min_gap"])
for q,e in E.items(): e["left"]=float(rms(R-e["w"])/FLOOR)
for q in OLD: INFO[q]["old_final"]=OLD[q][max(OLD[q])]; INFO[q]["old_best"]=min(OLD[q].values())
json.dump(INFO,open(f"{OUT}/results.json","w"),indent=1)
# ---------- Fig 7 style curves
fig,ax=plt.subplots(figsize=(3.6,5.4))
for q in STY:
    N=sorted(OLD[q]); ax.plot(N,[OLD[q][n] for n in N],color=STY[q],lw=1.0,ls=(0,(4,2)),alpha=0.9)
    N=sorted(NEW[q]); ax.plot(N,[NEW[q][n] for n in N],color=STY[q],lw=1.5,label=NAME[q])
ax.axhline(1,color="green",lw=1.5); ax.axhline(237.3,color="0.5",ls=":",lw=0.8); ax.axhline(88.8,color="0.6",ls=":",lw=0.8)
ax.text(1.15,1.12,"CMB flat, 1×",color="green",fontsize=6.5,va="bottom"); ax.text(1.15,250,"start 237×",color="0.4",fontsize=6.5,va="bottom"); ax.text(1.15,93,"four clumps 89×",color="0.45",fontsize=6.5,va="bottom")
ax.set_xscale("log"); ax.set_yscale("log"); xmax=max(max(NEW[q]) for q in NEW); ax.set_xlim(1,xmax*1.6)
ymin=min(0.5,0.5*min(min(NEW[q].values()) for q in NEW)); ymax=max(2000,1.5*max(max(NEW[q].values()) for q in NEW)); ax.set_ylim(ymin,ymax)
ax.set_xlabel("bodies placed"); ax.set_ylabel(r"leftover ($\times$ floor)"); ax.set_title("Sequential placement, core and our sphere open",fontsize=8)
hd=[Line2D([],[],color=STY[q],lw=1.5,label=NAME[q]) for q in STY]+[Line2D([],[],color="0.3",lw=1.5,label="open region (this run)"),Line2D([],[],color="0.3",lw=1.0,ls=(0,(4,2)),label="old runs (14.3–42.9 Gpc, core blocked)")]
ax.legend(handles=hd,loc="upper center",bbox_to_anchor=(0.5,-0.13),frameon=False,ncol=1,fontsize=6.5)
fig.savefig(f"{OUT}/fig7_open_curves.pdf",bbox_inches="tight"); fig.savefig(f"{OUT}/fig7_open_curves.png",bbox_inches="tight",dpi=200); plt.close(fig)
# ---------- Fig 8 style universes
SLAB=1.5; order=sorted(E,key=lambda q:-E[q]["left"])
fig=plt.figure(figsize=(7.1,10.6)); GO=fig.add_gridspec(3,2,hspace=0.30,wspace=0.25,left=0.07,right=0.98,top=0.965,bottom=0.09)
X0,X1,Y0,Y1=-106,48,-64,48; Ld=np.arange(360)
for i,q in enumerate(order):
    e=E[q]; col=STY[q]; sub=GO[i//2,i%2].subgridspec(2,1,height_ratios=[1.7,0.95],hspace=0.62); ax=fig.add_subplot(sub[0])
    ax.add_patch(Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="#e6b0aa",lw=0.6,ls="--",zorder=0))
    ax.text(C_G[0],C_G[1]-A_IN_G+4,"CoE core (open)",color="#a93226",fontsize=6,ha="center",va="bottom",zorder=6,bbox=dict(fc="white",ec="none",alpha=0.8,pad=0.5))
    ax.add_patch(Circle((0,0),S_MAX,fill=False,ec="0.6",lw=0.6,ls=":",zorder=1))
    Xv,rv,Xc=e["Xv"],e["rv"],e["Xc"]
    if e["slab"]: sv=np.abs(Xv[:,2])<SLAB; sc=np.abs(Xc[:,2])<SLAB; Xv,rv,Xc=Xv[sv],rv[sv],Xc[sc]
    ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(Xv,rv)],fc="#aed6f1",ec="#2471a3",lw=0.25,zorder=2))
    if len(Xc): ax.scatter(Xc[:,0],Xc[:,1],s=1.2 if e["slab"] else 6,c="#c0392b",lw=0,zorder=3)
    ax.add_patch(Circle((0,0),RMIN,fill=False,ec="#1f4e79",lw=1.4,zorder=4)); ax.plot(0,0,"k+",ms=5,mew=1,zorder=5)
    ax.set_xlim(X0,X1); ax.set_ylim(Y0,Y1); ax.set_aspect("equal"); ax.tick_params(labelsize=6)
    ax.set_xlabel("x (Gpc)",fontsize=6.5,labelpad=1); ax.set_ylabel("y (Gpc)",fontsize=6.5,labelpad=1)
    for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(2.2)
    tcol="#b48fc9" if q=="seqB_shallow" else col
    ax.set_title(f"{i+1}. {NAME[q]}\n{e['n']}"+("  (slab |z|<1.5 Gpc)" if e["slab"] else ""),color=tcol,fontsize=7.2,fontweight="bold")
    r_=fig.add_subplot(sub[1]); res=(R-e["w"])/FLOOR
    r_.axhspan(-1,1,color="green",alpha=0.25,lw=0); r_.axhline(1,color="green",lw=0.7); r_.axhline(-1,color="green",lw=0.7); r_.axhline(0,color="0.6",lw=0.5); r_.plot(Ld,res,color=col,lw=1.1)
    lim=1.15*np.abs(res).max(); r_.set_ylim(-lim,lim); r_.set_xlim(0,360); r_.set_xticks([0,90,180,270,360]); r_.tick_params(labelsize=6)
    r_.set_title(f"leftover {e['left']:,.1f}× the CMB level",fontsize=7.2,color=tcol,fontweight="bold")
    r_.set_xlabel("longitude ℓ (°) around the ring",fontsize=6.3,labelpad=1); r_.set_ylabel("× $10^{-5}$",fontsize=6.5,labelpad=1)
    for s_ in r_.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(1.2)
hd=[Line2D([],[],marker="o",ls="",color="#c0392b",ms=3.5,label="cluster"),Line2D([],[],marker="o",ls="",mfc="#aed6f1",mec="#2471a3",ms=7,label="void (to scale)"),
    Line2D([],[],color="#1f4e79",lw=1.4,label="CMB ring / our sphere (14.3 Gpc), open inside"),plt.Rectangle((0,0),1,1,fc="#fadbd8",ec="#e6b0aa",ls="--",label="CoE core (open to bodies)"),
    Line2D([],[],color="0.6",lw=0.6,ls=":",label="reach outside the core (42.9 Gpc)"),plt.Rectangle((0,0),1,1,fc="green",alpha=0.25,ec="none",label="±1× (CMB level)")]
fig.legend(handles=hd,loc="lower center",ncol=3,frameon=False,fontsize=6.6,bbox_to_anchor=(0.5,0.0))
fig.savefig(f"{OUT}/fig8_open_universes.pdf",bbox_inches="tight"); fig.savefig(f"{OUT}/fig8_open_universes.png",bbox_inches="tight",dpi=200); plt.close(fig)
for q in order: print(q,round(E[q]["left"],2),INFO[q])

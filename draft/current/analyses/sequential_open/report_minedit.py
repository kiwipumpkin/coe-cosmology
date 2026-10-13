"""Report for the minimal-edit test (minedit.py): curves, move/resize stats, where the edited bodies are, side views."""
import os,json,glob,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle,FancyArrowPatch
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
import pickle
from common import *
plt.rcParams.update({"font.size":7.5})
RUNS=[("b_random_s1_tidalcore","random s1","#95a5a6"),("b_random_s2_tidalcore","random s2","#7f8c8d"),("b_random_s3_tidalcore","random s3","#566573"),
      ("b_foam_rim_s1_tidalcore","foam s1","#85c1e9"),("b_foam_rim_s2_tidalcore","foam s2","#5dade2"),("b_foam_rim_s3_tidalcore","foam s3","#2e86c1")]
RUNS=[(t+sf,nm+lab,c) for sf,lab in (("_nr"," (near-ring bodies movable)"),(""," (near-ring bodies frozen)")) for t,nm,c in RUNS]
BASE=lambda t:t.replace("_nr","")
def dring(X): X=np.atleast_2d(X); return np.sqrt((np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2)
D={}
for t,nm,c in RUNS:
    fj=f"{OUT}/minedit/{t}.json"
    if os.path.exists(fj): D[t]=json.load(open(fj)); D[t]["final"]=True
    elif os.path.exists(f"{OUT}/ckpt/minedit_{t}.pkl"): D[t]=dict(curve=pickle.load(open(f"{OUT}/ckpt/minedit_{t}.pkl","rb"))["curve"],final=False)
runs=[r for r in RUNS if r[0] in D]
fig,ax=plt.subplots(figsize=(7,7.5))
for t,nm,c in runs:
    cv=np.array(D[t]["curve"]); x=np.maximum(cv[:,0],0.8); ax.plot(x,cv[:,1],"-" if t.endswith("_nr") else "--",color=c,lw=1.4,label=nm+("" if D[t]["final"] else " (running)"))
ax.axhline(1,color="green",lw=1.5); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("number of existing bodies moved / resized (log; 0 plotted at 0.8)"); ax.set_ylabel("leftover (× 1×)")
ax.set_title("Minimal edit: fewest existing bodies changed to flatten the ring (no new bodies)",fontsize=8.5)
fig.legend(handles=[l for l in ax.get_lines() if not l.get_label().startswith('_')]+[Line2D([],[],color="green",lw=1.5,label="1× (CMB level)")],loc="lower center",ncol=2,frameon=False,fontsize=7); fig.tight_layout(rect=(0,0.16,1,1))
fig.savefig(f"{OUT}/minedit_curve.png",dpi=200,bbox_inches="tight"); plt.close(fig)
fin=[r for r in runs if D[r[0]]["final"]]
if fin:
    fig,axs=plt.subplots(2,2,figsize=(7.2,6.4))
    for t,nm,c in fin:
        d=D[t]; ch=np.array(d["changed"]); mv=np.array(d["move"])[ch]; iv=np.array(d["is_void"])[ch]
        ro=np.array(d["r_old"])[ch]; rn=np.array(d["r_new"])[ch]; Mo=np.array(d["M_old"])[ch]; Mn=np.array(d["M_new"])[ch]
        axs[0,0].hist(mv,bins=np.linspace(0,15,31),histtype="step",color=c,lw=1.2,ls="-" if t.endswith("_nr") else "--")
        ratio=np.where(iv,rn/ro,Mn/Mo); axs[0,1].hist(np.log10(ratio),bins=np.linspace(-1.05,1.05,22),histtype="step",color=c,lw=1.2,ls="-" if t.endswith("_nr") else "--")
        axs[1,0].hist(np.array(d["dist_us_new"])[ch],bins=np.linspace(0,30,31),histtype="step",color=c,lw=1.2,ls="-" if t.endswith("_nr") else "--")
        axs[1,1].hist(np.array(d["dist_ring_new"])[ch],bins=np.linspace(0,25,26),histtype="step",color=c,lw=1.2,ls="-" if t.endswith("_nr") else "--")
    f0=np.load(f"{OUT}/final_{BASE(fin[0][0])}.npz"); Xa=f0["Xp"]
    axs[1,0].hist(np.linalg.norm(Xa,axis=1),bins=np.linspace(0,30,31),histtype="stepfilled",color="green",alpha=0.2,weights=np.full(len(Xa),np.mean([np.sum(D[t]["changed"]) for t,_,_ in fin])/len(Xa)))
    axs[1,1].hist(dring(Xa),bins=np.linspace(0,25,26),histtype="stepfilled",color="green",alpha=0.2,weights=np.full(len(Xa),np.mean([np.sum(D[t]["changed"]) for t,_,_ in fin])/len(Xa)))
    axs[1,1].axvline(3,color="red",ls=":",lw=0.8)
    axs[0,0].set_xlabel("how far each edited body moved (Gpc)"); axs[0,1].set_xlabel("log10 size change (void radius ratio or cluster mass ratio)")
    axs[1,0].set_xlabel("edited body's distance from us (Gpc, new position)"); axs[1,1].set_xlabel("edited body's distance from the ring (Gpc, new)")
    for a in axs.ravel(): a.set_ylabel("bodies")
    hd=[Line2D([],[],color=c,ls="-" if t.endswith("_nr") else "--",label=nm) for t,nm,c in fin]+[plt.Rectangle((0,0),1,1,fc="green",alpha=0.2,label="all web bodies (scaled to the same count)")]
    fig.legend(handles=hd,loc="upper center",ncol=2,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.0)); fig.tight_layout(); fig.savefig(f"{OUT}/minedit_changes.png",dpi=200,bbox_inches="tight"); plt.close(fig)
    # side views: up to 4 finished runs
    sv=sorted(fin,key=lambda r:(not r[0].endswith("_nr"),D[r[0]]["left"]))[:4]; fig=plt.figure(figsize=(10.5,5.6*int(np.ceil(len(sv)/2)))); G=fig.add_gridspec(int(np.ceil(len(sv)/2)),2,hspace=0.5,wspace=0.3)
    for q,(t,nm,c) in enumerate(sv):
        d=D[t]; f=np.load(f"{OUT}/final_{BASE(t)}.npz"); v=f["kp"]=="void"; Xp=f["Xp"]; rp=f["rp"]; sub=G[q//2,q%2].subgridspec(2,1,height_ratios=[1.7,0.8],hspace=0.45); ax=fig.add_subplot(sub[0])
        cp=Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="none",zorder=0); ax.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=ax.transData))
        ed=np.zeros(len(Xp),bool); ed[np.array(d["idx"])[np.array(d["changed"])]]=True
        sl=(np.abs(Xp[:,2])<1.5)&~ed
        ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(Xp[sl&v],rp[sl&v])],fc="#eaf2f8",ec="#d4e6f1",lw=0.15,zorder=1))
        ax.scatter(Xp[sl&~v,0],Xp[sl&~v,1],s=0.4,c="#f2d7d5",lw=0,zorder=1)
        ch=np.array(d["changed"]); Xo=np.array(d["X_old"])[ch]; Xn=np.array(d["X_new"])[ch]; iv=np.array(d["is_void"])[ch]; rn=np.array(d["r_new"])[ch]; ro=np.array(d["r_old"])[ch]
        ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(Xo[iv],ro[iv])],fc="none",ec="0.55",lw=0.6,ls="--",zorder=3))
        ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(Xn[iv],rn[iv])],fc="#2471a3",alpha=0.35,ec="#1a5276",lw=0.8,zorder=4))
        ax.scatter(Xo[~iv,0],Xo[~iv,1],s=6,facecolor="none",edgecolor="0.4",lw=0.5,zorder=3); ax.scatter(Xn[~iv,0],Xn[~iv,1],s=10,c="#c0392b",edgecolor="k",lw=0.3,zorder=5)
        for a_,b_ in zip(Xo,Xn):
            if np.linalg.norm(b_-a_)>0.2: ax.annotate("",xy=b_[:2],xytext=a_[:2],arrowprops=dict(arrowstyle="-|>",color="k",lw=0.5,mutation_scale=5),zorder=6)
        th=np.linspace(0,2*np.pi,300)
        for rr in (RMIN-3,RMIN+3): ax.plot(rr*np.cos(th),rr*np.sin(th),color="red",lw=0.5,ls=":",zorder=2)
        ax.add_patch(Circle((0,0),RMIN,fill=False,ec="#1f4e79",lw=1.3,zorder=5)); ax.plot(0,0,"k+"); frame30(ax); ax.tick_params(labelsize=5.5)
        for s_ in ax.spines.values(): s_.set_edgecolor(c); s_.set_linewidth(2)
        ax.set_title(f"{nm}:\n{d['n_edited']} of {d['n_bodies']:,} bodies changed ({100*d['frac_edited']:.2f}%), all z projected\nleftover {d['before']['total']:.0f}× → {d['left']:.2f}×; median move {np.median(np.array(d['move'])[ch]):.1f} Gpc",fontsize=6.5,color=c,fontweight="bold")
        r_=fig.add_subplot(sub[1]); r_.axhspan(-1,1,color="green",alpha=0.3,lw=0); r_.axhline(0,color="0.6",lw=0.5)
        f0=np.load(f"{OUT}/final_{BASE(t)}.npz"); before=(R-f0["w"])/FLOOR; after=np.array(d["resid"])/FLOOR
        r_.plot(np.arange(360),before,color="0.6",lw=0.8); r_.plot(np.arange(360),after,color=c,lw=1.0)
        lim=1.3*max(np.abs(before).max(),2); r_.set_yscale("symlog",linthresh=1,linscale=0.6); r_.set_ylim(-lim,lim); r_.set_xlim(0,360); r_.set_xticks([0,90,180,270,360]); r_.tick_params(labelsize=5.5)
        fs=d["after"]; r_.set_title(f"before (grey) and after: m1–12 {np.hypot(fs['rms_m1_5'],fs['rms_m6_12']):.2f}, 13+ {fs['rms_13up']:.2f}, symmetric {fs['sym']*100:.0f}%",fontsize=6.3,color=c)
        r_.set_xlabel("ℓ (°) around the ring",fontsize=5.8); r_.set_ylabel("× $10^{-5}$ (symlog)",fontsize=5.5)
    hd=[Line2D([],[],marker="o",ls="",mfc="none",mec="0.4",label="changed body, old position/size"),Line2D([],[],marker="o",ls="",mfc="#c0392b",mec="k",label="changed cluster, new"),
        plt.Rectangle((0,0),1,1,fc="#2471a3",alpha=0.35,ec="#1a5276",label="changed void, new (to scale)"),Line2D([],[],color="k",lw=0.6,marker=">",label="move (old → new)"),
        Line2D([],[],marker="o",ls="",mfc="#eaf2f8",mec="#d4e6f1",label="unchanged web (slab |z|<1.5)"),Line2D([],[],color="red",ls=":",lw=0.6,label="±3 Gpc from the ring"),
        Line2D([],[],color="#1f4e79",lw=1.3,label="CMB ring / our sphere"),plt.Rectangle((0,0),1,1,fc="#fadbd8",label="core cap inside 30 Gpc"),Line2D([],[],color="#a93226",marker=">",label="direction to the CoE"),
        Line2D([],[],color="0.6",lw=0.8,label="ring residual before"),plt.Rectangle((0,0),1,1,fc="green",alpha=0.3,label="±1× (CMB level)")]
    fig.legend(handles=hd,loc="upper center",ncol=4,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.0)); fig.savefig(f"{OUT}/minedit_side_views.png",dpi=170,bbox_inches="tight"); plt.close(fig)
S={t:(dict(n_bodies=D[t]["n_bodies"],edited=D[t]["n_edited"],frac=round(D[t]["frac_edited"],4),before=round(D[t]["before"]["total"],1),left=round(D[t]["left"],2),reached=D[t]["reached"],
          m1_12=round(float(np.hypot(D[t]["after"]["rms_m1_5"],D[t]["after"]["rms_m6_12"])),2),m13=round(D[t]["after"]["rms_13up"],2),sym=round(D[t]["after"]["sym"],2),
          median_move=round(float(np.median(np.array(D[t]["move"])[np.array(D[t]["changed"])])),2),net_mass_change=D[t]["net_mass_change"]) if D[t]["final"] else dict(running=True,curve_last=D[t]["curve"][-1])) for t,_,_ in runs}
json.dump(S,open(f"{OUT}/minedit/summary.json","w"),indent=1)
for k,v in S.items(): print(k,v)

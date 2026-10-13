"""Default edit figure (Chris, 10 Oct): plain before/after pairs. For each run: 'before' side view (web as grown) and 'after' side view
(edited web) side by side on the same axes (modeled region only, frame30), no arrows / deletion markers; ring residual under each
(green ±1× band, symlog); legend below."""
import json,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from common import *
def layout(base,d=None):
    f=np.load(f"{OUT}/final_{base}.npz"); X=f["Xp"].astype(float).copy(); r=f["rp"].astype(float).copy(); v=(f["kp"]=="void"); M=f["Mp"].astype(float).copy()
    if d is not None:
        i=np.array(d["idx"],int); X[i]=np.array(d["X_new"]); r[i]=np.array(d["r_new"]); M[i]=np.array(d["M_new"])
    return X,r,v,M,f
HL="#1e8449"; HLE="#0b5345"; DEL="#8e44ad"; BAND="#d5f5e3"   # moved = dark green; ±1x band = pale mint
def draw_map(ax,X,r,v,M,title,col,hl=None,dl=None):
    """hl: bodies that move -> green (before: original spot, after: new spot); resized-only bodies keep red/blue; dl: deleted (BEFORE only) -> purple"""
    n=len(X); hl=np.zeros(n,bool) if hl is None else hl; dl=np.zeros(n,bool) if dl is None else dl
    cp=Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="none",zorder=0); ax.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=ax.transData))
    al=r>0; sl=(np.abs(X[:,2])<1.5)&al&(np.linalg.norm(X,axis=1)<=30); pl=al&(np.linalg.norm(X,axis=1)<=30)   # unchanged: slab; highlighted: all z
    u=sl&~hl&~dl
    ax.add_collection(PatchCollection([Circle(x[:2],rr) for x,rr in zip(X[u&v],r[u&v])],fc="#aed6f1",ec="#5dade2",lw=0.2,alpha=0.6,zorder=1))
    c=u&~v; ax.scatter(X[c,0],X[c,1],s=np.clip(6*(M[c]/1.4e50)**(2/3),0.5,40),c="#c0392b",lw=0,alpha=0.8,zorder=2)
    for m_,fc,ec,z in ((pl&hl,HL,HLE,4),(pl&dl,DEL,"#4a235a",4)):
        ax.add_collection(PatchCollection([Circle(x[:2],rr) for x,rr in zip(X[m_&v],r[m_&v])],fc=fc,ec=ec,lw=0.5,alpha=0.55,zorder=z))
        cc=m_&~v; ax.scatter(X[cc,0],X[cc,1],s=np.clip(9*(np.abs(M[cc])/1.4e50)**(2/3),2,50),c=fc,edgecolor=ec,lw=0.3,zorder=z+1)
    th=np.linspace(0,2*np.pi,300)
    for rr in (RMIN-3,RMIN+3): ax.plot(rr*np.cos(th),rr*np.sin(th),color="red",lw=0.5,ls=":",zorder=3)
    ax.add_patch(Circle((0,0),RMIN,fill=False,ec="#1f4e79",lw=1.3,zorder=6)); ax.plot(0,0,"k+",zorder=7); frame30(ax); ax.tick_params(labelsize=5.5)
    for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(1.6)
    ax.set_title(title,fontsize=6.5,color=col,fontweight="bold")
def draw_res(ax,res,col,lim,title):
    ax.axhspan(-1,1,color=BAND,alpha=1.0,lw=0,zorder=0); ax.axhline(1,color="#58d68d",lw=0.5,ls="--"); ax.axhline(-1,color="#58d68d",lw=0.5,ls="--"); ax.axhline(0,color="0.6",lw=0.5); ax.plot(np.arange(360),res,color=col,lw=0.9)
    ax.set_yscale("symlog",linthresh=1,linscale=0.6); ax.set_ylim(-lim,lim); ax.set_xlim(0,360); ax.set_xticks([0,90,180,270,360]); ax.tick_params(labelsize=5.5)
    ax.set_xlabel("ℓ (°) around the ring",fontsize=5.8); ax.set_ylabel("× $10^{-5}$ (symlog)",fontsize=5.5); ax.set_title(title,fontsize=6.2,color=col)
def figure(items,path):
    """items: list of (json_path, label, colour)"""
    anydel=any(any(json.load(open(js)).get("deleted",[])) for js,_,_ in items)
    nr=len(items); fig=plt.figure(figsize=(10.5,6.6*nr)); G=fig.add_gridspec(nr,1,hspace=0.32)
    for k,(js,lab,col) in enumerate(items):
        d=json.load(open(js)); base=d["tag"].split("_clr")[0].split("_nodel")[0].replace("_pruned","").replace("_nr","")
        Xb,rb,vb,Mb,f=layout(base); Xa,ra,va,Ma,_=layout(base,d)
        sub=G[k].subgridspec(2,2,height_ratios=[1.7,0.75],hspace=0.42,wspace=0.22)
        before=(R-f["w"])/FLOOR; after=np.array(d["resid"])/FLOOR; lim=1.3*max(np.abs(before).max(),2); sb=d["before"]; sa=d["after"]
        idx=np.array(d["idx"],int); de=np.array(d["deleted"],bool); ch=np.array(d["moved"],bool)
        HLm=np.zeros(len(Xb),bool); HLm[idx[ch&~de]]=True; DLm=np.zeros(len(Xb),bool); DLm[idx[de]]=True
        nd=d.get("n_deleted",0); nm=d.get("n_moved",0); nz=d.get("n_resized",0)
        draw_map(fig.add_subplot(sub[0,0]),Xb,rb,vb,Mb,f"{lab}\nBEFORE: web as grown, {len(Xb):,} bodies\nleftover {sb['total']:.0f}×",col,hl=HLm,dl=DLm)
        draw_map(fig.add_subplot(sub[0,1]),Xa,ra,va,Ma,f"{lab}\nAFTER: {d['n_edited']} bodies edited ({100*d['frac_edited']:.1f}%): {nd} deleted, {nm} moved, {nz} resized\nleftover {d['left']:.2f}×",col,hl=HLm)
        f2=lambda s_: f"m1–12 {np.hypot(s_['rms_m1_5'],s_['rms_m6_12']):.2f}×, 13+ {s_['rms_13up']:.2f}×, symmetric {s_['sym']*100:.0f}%"
        draw_res(fig.add_subplot(sub[1,0]),before,"0.35",lim,"before: "+f2(sb)); draw_res(fig.add_subplot(sub[1,1]),after,col,lim,"after: "+f2(sa))
    hd=[plt.Rectangle((0,0),1,1,fc="#aed6f1",ec="#5dade2",alpha=0.6,label="void not moved (to scale, slab |z|<1.5 Gpc; resized-only included)"),Line2D([],[],marker="o",ls="",color="#c0392b",label="cluster not moved (size ∝ mass^2/3; resized-only included)"),
        plt.Rectangle((0,0),1,1,fc=HL,ec=HLE,alpha=0.55,label="void that moves (before: original spot & size; after: new spot & size; all z)"),Line2D([],[],marker="o",ls="",mfc=HL,mec=HLE,label="cluster that moves (before: original spot; after: new spot; all z)"),
        plt.Rectangle((0,0),1,1,fc=DEL,ec="#4a235a",alpha=0.55,label="void that is deleted (before view only)"),Line2D([],[],marker="o",ls="",mfc=DEL,mec="#4a235a",label="cluster that is deleted (before view only)"),
        Line2D([],[],color="#1f4e79",lw=1.3,label="CMB ring / our sphere"),Line2D([],[],color="red",ls=":",lw=0.6,label="±3 Gpc from the ring"),
        plt.Rectangle((0,0),1,1,fc="#fadbd8",label="core cap inside 30 Gpc"),Line2D([],[],color="0.4",ls=":",label="30 Gpc (modeled region)"),
        Line2D([],[],color="#a93226",marker=">",label="direction to the CoE"),plt.Rectangle((0,0),1,1,fc=BAND,ec="#58d68d",ls="--",lw=0.5,label="±1× CMB band (pale, dashed edges)")]
    if not anydel: hd=[h for h in hd if "deleted" not in h.get_label()]
    fig.legend(handles=hd,loc="upper center",ncol=3,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.0)); fig.savefig(path,dpi=160,bbox_inches="tight"); plt.close(fig)

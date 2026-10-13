"""Default edit figure (Chris, 10 Oct): ONE PNG PER RUN (run_figure): before | after side views, residual below, legend below.
(figure() = older multi-run stack, kept for reference.)
One run per row (~4:3): BEFORE map | AFTER map on top, ring residual (before grey, after coloured) full width below; legend below the figure.
Colours: clusters red, voids blue; bodies that MOVE are dark green (BEFORE: at their original spot, AFTER: at their new spot);
resized-only bodies stay red/blue; deleted bodies purple in BEFORE only; ±1x CMB band pale mint with dashed edges. No arrows.
Fonts >= ~18 pt on a 12-inch-wide figure (base FS)."""
import json,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from common import *
FS=12; FT=13; FK=11; FL=11          # moderate fonts (Chris, 10 Oct): one run per image instead of big fonts
HL="#1e8449"; HLE="#0b5345"; DEL="#8e44ad"; BAND="#d5f5e3"
def layout(base,d=None):
    f=np.load(f"{OUT}/final_{base}.npz"); X=f["Xp"].astype(float).copy(); r=f["rp"].astype(float).copy(); v=(f["kp"]=="void"); M=f["Mp"].astype(float).copy()
    if d is not None:
        i=np.array(d["idx"],int); X[i]=np.array(d["X_new"]); r[i]=np.array(d["r_new"]); M[i]=np.array(d["M_new"])
    return X,r,v,M,f
def draw_map(ax,X,r,v,M,title,col,hl=None,dl=None):
    n=len(X); hl=np.zeros(n,bool) if hl is None else hl; dl=np.zeros(n,bool) if dl is None else dl
    cp=Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="none",zorder=0); ax.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=ax.transData))
    al=r>0; sl=(np.abs(X[:,2])<1.5)&al&(np.linalg.norm(X,axis=1)<=30); pl=al&(np.linalg.norm(X,axis=1)<=30)
    u=sl&~hl&~dl
    ax.add_collection(PatchCollection([Circle(x[:2],rr) for x,rr in zip(X[u&v],r[u&v])],fc="#aed6f1",ec="#5dade2",lw=0.3,alpha=0.6,zorder=1))
    c=u&~v; ax.scatter(X[c,0],X[c,1],s=np.clip(10*(M[c]/1.4e50)**(2/3),1,60),c="#c0392b",lw=0,alpha=0.85,zorder=2)
    for m_,fc,ec,z in ((pl&hl,HL,HLE,4),(pl&dl,DEL,"#4a235a",4)):
        ax.add_collection(PatchCollection([Circle(x[:2],rr) for x,rr in zip(X[m_&v],r[m_&v])],fc=fc,ec=ec,lw=0.6,alpha=0.55,zorder=z))
        cc=m_&~v; ax.scatter(X[cc,0],X[cc,1],s=np.clip(16*(np.abs(M[cc])/1.4e50)**(2/3),4,80),c=fc,edgecolor=ec,lw=0.4,zorder=z+1)
    th=np.linspace(0,2*np.pi,300)
    for rr in (RMIN-3,RMIN+3): ax.plot(rr*np.cos(th),rr*np.sin(th),color="red",lw=0.9,ls=":",zorder=3)
    ax.add_patch(Circle((0,0),RMIN,fill=False,ec="#1f4e79",lw=2.2,zorder=6)); ax.plot(0,0,"k+",ms=12,mew=2,zorder=7); frame30(ax,lw=1.2)
    for t in ax.texts: t.set_fontsize(FK)
    ax.set_xticks([-30,0,30]); ax.set_yticks([-30,0,30]); ax.tick_params(labelsize=FK)
    for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(2.5)
    ax.set_title(title,fontsize=FT,color="k",linespacing=1.15)
def draw_res(ax,before,after,col,lim,tb,ta):
    ax.axhspan(-1,1,color=BAND,lw=0,zorder=0); ax.axhline(1,color="#58d68d",lw=1,ls="--"); ax.axhline(-1,color="#58d68d",lw=1,ls="--"); ax.axhline(0,color="0.6",lw=0.6)
    ax.plot(np.arange(360),before,color="0.55",lw=1.6,label=tb); ax.plot(np.arange(360),after,color="#1f3a93",lw=2.0,label=ta)
    ax.set_yscale("symlog",linthresh=1,linscale=0.6); ax.set_ylim(-lim,lim); ax.set_xlim(0,360); ax.set_xticks([0,90,180,270,360])
    ax.set_yticks([-1e3,-10,0,10,1e3]); ax.set_yticklabels(["−10³","−10","0","10","10³"]); ax.tick_params(labelsize=FK)
    ax.set_xlabel("position around the ring, ℓ (°)",fontsize=FS); ax.set_ylabel("residual (×10⁻⁵)",fontsize=FS)
    ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=1,frameon=False,fontsize=FK,handlelength=1.5,borderaxespad=0.2)
def figure(items,path):
    """items: list of (json_path, label, colour); one run per row"""
    anydel=any(any(json.load(open(js)).get("deleted",[])) for js,_,_ in items)
    nr=len(items); fig=plt.figure(figsize=(12,12.2*nr)); G=fig.add_gridspec(nr,1,hspace=0.16,bottom=0.012,top=0.995,left=0.09,right=0.98)
    for k,(js,lab,col) in enumerate(items):
        d=json.load(open(js)); base=d["tag"].split("_clr")[0].split("_nodel")[0].replace("_pruned","").replace("_nr","").replace("_entropy","")
        Xb,rb,vb,Mb,f=layout(base); Xa,ra,va,Ma,_=layout(base,d)
        sub=G[k].subgridspec(3,2,height_ratios=[0.06,1.0,0.42],hspace=0.62,wspace=0.12)
        before=(R-f["w"])/FLOOR; after=np.array(d["resid"])/FLOOR; lim=1.3*max(np.abs(before).max(),2); sb=d["before"]; sa=d["after"]
        idx=np.array(d["idx"],int); de=np.array(d["deleted"],bool); ch=np.array(d["moved"],bool)
        HLm=np.zeros(len(Xb),bool); HLm[idx[ch&~de]]=True; DLm=np.zeros(len(Xb),bool); DLm[idx[de]]=True
        nd=d.get("n_deleted",0); nm=d.get("n_moved",0); nz=d.get("n_resized",0)
        hx=fig.add_subplot(sub[0,:]); hx.axis("off"); hx.text(0.5,0.5,lab,ha="center",va="center",fontsize=FT+3,fontweight="bold",color=col if col not in ("#95a5a6","#717d7e","#85c1e9") else "#34495e")
        a0=fig.add_subplot(sub[1,0]); a1=fig.add_subplot(sub[1,1])
        draw_map(a0,Xb,rb,vb,Mb,f"BEFORE: as grown\n{len(Xb):,} bodies\nleftover {sb['total']:.0f}×",col,hl=HLm,dl=DLm)
        ed=f"{nm} moved, {nz} resized"+(f", {nd} deleted" if nd else "")
        draw_map(a1,Xa,ra,va,Ma,f"AFTER: {d['n_edited']} edited ({100*d['frac_edited']:.1f}%)\n{ed}\nleftover {d['left']:.2f}×",col,hl=HLm)
        f2=lambda s_: f"m1–12 {np.hypot(s_['rms_m1_5'],s_['rms_m6_12']):.2f}×, m13+ {s_['rms_13up']:.2f}×, sym {s_['sym']*100:.0f}%"
        draw_res(fig.add_subplot(sub[2,:]),before,after,col,lim,"before — "+f2(sb),"after — "+f2(sa))
    hd=[plt.Rectangle((0,0),1,1,fc="#aed6f1",ec="#5dade2",alpha=0.6,label="void, not moved (slab |z|<1.5)"),Line2D([],[],marker="o",ls="",ms=9,color="#c0392b",label="cluster, not moved"),
        plt.Rectangle((0,0),1,1,fc=HL,ec=HLE,alpha=0.55,label="void that moves (all z)"),Line2D([],[],marker="o",ls="",ms=9,mfc=HL,mec=HLE,label="cluster that moves (all z)"),
        plt.Rectangle((0,0),1,1,fc=DEL,ec="#4a235a",alpha=0.55,label="void deleted (BEFORE only)"),Line2D([],[],marker="o",ls="",ms=9,mfc=DEL,mec="#4a235a",label="cluster deleted (BEFORE only)"),
        Line2D([],[],color="#1f4e79",lw=2.2,label="CMB ring / our sphere"),Line2D([],[],color="red",ls=":",lw=1.2,label="±3 Gpc from the ring"),
        plt.Rectangle((0,0),1,1,fc="#fadbd8",label="CoE core cap (inside 30 Gpc)"),Line2D([],[],color="0.4",ls=":",lw=1.2,label="30 Gpc edge of model"),
        Line2D([],[],color="#a93226",marker=">",label="direction to the CoE"),plt.Rectangle((0,0),1,1,fc=BAND,ec="#58d68d",ls="--",lw=1,label="±1× CMB band (pale mint)")]
    if not anydel: hd=[h for h in hd if "deleted" not in h.get_label()]
    fig.legend(handles=hd,loc="upper center",ncol=2,frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.0),handlelength=1.6,columnspacing=1.5,
               title="moved bodies are green in both views; resized-only bodies keep red/blue",title_fontsize=FK-1)
    fig.savefig(path,dpi=110,bbox_inches="tight"); plt.close(fig)

def legend_handles(anydel):
    hd=[plt.Rectangle((0,0),1,1,fc="#aed6f1",ec="#5dade2",alpha=0.6,label="void, not moved (slab |z|<1.5 Gpc)"),Line2D([],[],marker="o",ls="",ms=6,color="#c0392b",label="cluster, not moved"),
        plt.Rectangle((0,0),1,1,fc=HL,ec=HLE,alpha=0.55,label="void that moves (all z)"),Line2D([],[],marker="o",ls="",ms=6,mfc=HL,mec=HLE,label="cluster that moves (all z)"),
        plt.Rectangle((0,0),1,1,fc=DEL,ec="#4a235a",alpha=0.55,label="void deleted (BEFORE only)"),Line2D([],[],marker="o",ls="",ms=6,mfc=DEL,mec="#4a235a",label="cluster deleted (BEFORE only)"),
        Line2D([],[],color="#1f4e79",lw=2,label="CMB ring / our sphere"),Line2D([],[],color="red",ls=":",lw=1,label="±3 Gpc from the ring"),
        plt.Rectangle((0,0),1,1,fc="#fadbd8",label="CoE core cap (inside 30 Gpc)"),Line2D([],[],color="0.4",ls=":",lw=1,label="30 Gpc edge of model"),
        Line2D([],[],color="#a93226",marker=">",label="direction to the CoE"),plt.Rectangle((0,0),1,1,fc=BAND,ec="#58d68d",ls="--",lw=1,label="±1× CMB band (pale mint)")]
    return hd if anydel else [h for h in hd if "deleted" not in h.get_label()]
def run_figure(js,lab,col,path):
    """one run: BEFORE | AFTER maps, residual around the ring below (before grey, after dark blue, pale-mint ±1x band), legend below"""
    d=json.load(open(js)); base=d["tag"].split("_clr")[0].split("_nodel")[0].replace("_pruned","").replace("_nr","").replace("_entropy","")
    Xb,rb,vb,Mb,f=layout(base); Xa,ra,va,Ma,_=layout(base,d)
    fig=plt.figure(figsize=(10,10)); G=fig.add_gridspec(2,2,height_ratios=[1.0,0.45],hspace=0.42,wspace=0.12,left=0.09,right=0.98,top=0.88,bottom=0.08)
    before=(R-f["w"])/FLOOR; after=np.array(d["resid"])/FLOOR; lim=1.3*max(np.abs(before).max(),2); sb=d["before"]; sa=d["after"]
    idx=np.array(d["idx"],int); de=np.array(d["deleted"],bool); ch=np.array(d["moved"],bool)
    HLm=np.zeros(len(Xb),bool); HLm[idx[ch&~de]]=True; DLm=np.zeros(len(Xb),bool); DLm[idx[de]]=True
    nd=d.get("n_deleted",0); nm=d.get("n_moved",0); nz=d.get("n_resized",0); ed=f"{nm} moved, {nz} resized"+(f", {nd} deleted" if nd else "")
    fig.suptitle(lab,fontsize=FT+2,fontweight="bold",color=col if col not in ("#95a5a6","#717d7e","#85c1e9") else "#34495e",y=0.975)
    draw_map(fig.add_subplot(G[0,0]),Xb,rb,vb,Mb,f"BEFORE: as grown, {len(Xb):,} bodies\nleftover {sb['total']:.0f}×",col,hl=HLm,dl=DLm)
    draw_map(fig.add_subplot(G[0,1]),Xa,ra,va,Ma,f"AFTER: {d['n_edited']} edited ({100*d['frac_edited']:.1f}%), {ed}\nleftover {d['left']:.2f}×",col,hl=HLm)
    f2=lambda s_: f"m1–12 {np.hypot(s_['rms_m1_5'],s_['rms_m6_12']):.2f}×, m13+ {s_['rms_13up']:.2f}×, symmetric {s_['sym']*100:.0f}%"
    ax=fig.add_subplot(G[1,:]); draw_res(ax,before,after,col,lim,"before — "+f2(sb),"after — "+f2(sa))
    ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=2,frameon=False,fontsize=FK,handlelength=1.5,borderaxespad=0.2)
    fig.legend(handles=legend_handles(nd>0),loc="upper center",ncol=3,frameon=False,fontsize=FL,bbox_to_anchor=(0.5,0.0),handlelength=1.6,columnspacing=1.2,
               title="moved bodies are green in both views; resized-only bodies keep red/blue",title_fontsize=FK)
    fig.savefig(path,dpi=130,bbox_inches="tight"); plt.close(fig)

def residual_below(fig,slot,d,f,col="#1f3a93"):
    """DEFAULT for every before/after side-view figure (Chris, 10 Oct): residual-around-the-ring panel below the maps.
    slot = a gridspec cell; d = run json; f = loaded web npz (layout()[4])."""
    before=(R-f["w"])/FLOOR; after=np.array(d["resid"])/FLOOR; lim=1.3*max(np.abs(before).max(),2); sb=d["before"]; sa=d["after"]
    f2=lambda s_,L: f"{L:.4g}× — m1–12 {np.hypot(s_['rms_m1_5'],s_['rms_m6_12']):.2f}×, m13+ {s_['rms_13up']:.2f}×, symmetric {s_['sym']*100:.0f}%"
    ax=fig.add_subplot(slot); draw_res(ax,before,after,col,lim,"before "+f2(sb,sb["total"]),"after "+f2(sa,d["left"]))
    ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=2,frameon=False,fontsize=FK,handlelength=1.5,borderaxespad=0.2); return ax
def change_map(ax,GX,GY,D,lim,linthresh,hatch=False):
    """DEFAULT area-change map (Chris, 10 Oct, final): diverging blue-white-red (RdBu_r): BLUE = shrank (compressed, 'blueshift'),
    WHITE = no change, RED = grew (expanded, 'redshift'). No hatch; points with no face in one view (non-crossing mesh only, or the
    rim outside a triangulation's hull) show a plain light-grey background."""
    from matplotlib.colors import SymLogNorm
    ax.set_facecolor("#e6e6e6")
    if hatch: ax.add_patch(plt.Rectangle((-40,-40),80,80,fc="#d9d9d9",ec="#9a9a9a",hatch="///",lw=0,zorder=-1))
    return ax.pcolormesh(GX,GY,D,cmap="RdBu_r",norm=SymLogNorm(linthresh=linthresh,vmin=-lim,vmax=lim),shading="auto",zorder=0)
CHANGE_LABEL="blue = compressed (smaller), red = expanded (bigger), white = no change"
NOFACE=dict(fc="#e6e6e6",ec="0.6")

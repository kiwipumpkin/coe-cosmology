"""Side views (ring plane seen face-on) of chosen layouts with network links, our sphere, and the leftover around the ring."""
import json,numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from jn_core import *
plt.rcParams.update({"font.size":12,"axes.titlesize":12.5})
def mix(k): return f"{k.count('C')} masses + {k.count('V')} voids"
def side_view(ax,r,lim=None):
    k=np.array(list(r["kinds"])); X=np.array(r["X"]); M=np.array(r["M"]); Rv=np.array(r["Rv"]); t=np.linspace(0,2*np.pi,400)
    far=np.abs(X[:,:2]).max()+4
    lim=lim or float(min(70,max(24,far)))
    ax.set_aspect("equal"); ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
    ax.plot(C_G[0]+A_IN_G*np.cos(t),C_G[1]+A_IN_G*np.sin(t),color="0.6",ls="--",lw=1)
    ax.plot(RV*np.cos(t),RV*np.sin(t),"k",lw=1.6)
    for a,b in r["props"]["links"]: ax.plot([X[a,0],X[b,0]],[X[a,1],X[b,1]],color="0.35",lw=0.9,alpha=0.7,zorder=2)
    for i in range(len(k)):
        if k[i]=="C": ax.scatter(X[i,0],X[i,1],s=25+220*abs(M[i])/2e52,c="#c0392b",edgecolors="k",lw=0.5,zorder=4)
        else: ax.add_patch(plt.Circle(X[i,:2],Rv[i],color="#2471a3",alpha=0.30+0.5*abs(r["delta"][i])/0.7,zorder=1))
    off=np.abs(X[:,:2]).max(1)>lim
    if off.any(): ax.text(0.02,0.02,f"{off.sum()} bodies beyond the frame",transform=ax.transAxes,fontsize=10,color="0.3")
    ax.plot(0,0,"k+",ms=13,mew=2); ax.set_xlabel("Gpc"); ax.tick_params(labelsize=10)
def ring(ax,r):
    deg=np.degrees(TH); v=np.array(r["resid"])/FLOOR
    ax.plot(deg,v,color="#34495e",lw=1.6); ax.axhline(0,color="0.5",lw=0.8)
    ax.axhspan(-1,1,color="green",alpha=0.22,lw=0,label="±1× CMB band")
    ax.set_xlim(0,360); ax.set_xticks([0,90,180,270,360]); ax.set_xlabel("position around the ring (deg)"); ax.set_ylabel("leftover (× floor)")
    ax.tick_params(labelsize=10)
def sheet(items,title,fn,ncol=2):
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    n=len(items); nrow=(n+ncol-1)//ncol
    fig=plt.figure(figsize=(7.4*ncol,10.2*nrow)); H=10.2*nrow; outer=GridSpec(nrow,ncol,hspace=0.36,wspace=0.25,top=1-1.9/H,bottom=1.3/H)
    for i,(lab,r) in enumerate(items):
        row,col=divmod(i,ncol); g=GridSpecFromSubplotSpec(2,1,subplot_spec=outer[row,col],height_ratios=[3.4,1],hspace=0.28)
        a=fig.add_subplot(g[0]); b=fig.add_subplot(g[1])
        import textwrap
        t1="\n".join(textwrap.wrap(lab.split("\n")[0],52)); t2=lab.split("\n")[1] if "\n" in lab else ""
        side_view(a,r); a.set_title(f"{t1}\n{t2}\n{mix(r['kinds'])}, leftover {r['left']:.2f}×\n{r['props']['n_inside']} inside our sphere, "
                                    f"{len(r['kinds'])-r['props']['n_inside']} outside ({r['props']['n_core']} in the CoE core)",fontsize=11)
        ring(b,r)
    fig.suptitle(title,fontsize=15,y=1-0.25/H)
    fig.text(0.5,0.15/H,"Side view = the ring plane seen face-on (bodies above/below it are projected). + = us, black circle = our 14.3 Gpc sphere,\n"
             "grey dashed = CoE core edge, red dots = masses (size ~ mass), blue discs = voids (darker = deeper), grey lines = network links.\n"
             "Strip under each view: leftover around the ring on its own scale; green band = ±1× (CMB flat).",ha="center",fontsize=11.5)
    fig.savefig(fn,dpi=85); plt.close(fig)

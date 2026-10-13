"""net_force_arrows.png (+ _zoom): each slab body (|z|<1.5 Gpc) replaced by an arrow centred on it, along its in-plane net force
(Delaunay-neighbour pair forces from pair_energy.py). Length: voids = body diameter 3 Gpc, clusters = 1.0 Gpc (their 0.32 Gpc diameter is too
small to see), each times (0.6 + 0.4 q), q = rank of log|F| within its type (0..1): a mild log-|F| scaling. Arrows are the explicit exception
to the no-arrows rule (force maps only)."""
import numpy as np
src=open("pair_energy.py").read(); exec(src[:src.index("# ---------- (1)")])
COE="#641e16"
def coe(ax,lim=None):
    u=C_G[:2]/np.linalg.norm(C_G[:2]); L=30.0
    ax.annotate("",xy=(L*u[0],L*u[1]),xytext=(0,0),arrowprops=dict(arrowstyle="-|>,head_length=1.0,head_width=0.5",color=COE,lw=3.6,shrinkA=0,shrinkB=0),zorder=12)
    p=19*u+np.array([-u[1],u[0]])*2.8; ax.text(p[0],p[1],"to CoE",color=COE,fontsize=11,fontweight="bold",ha="center",va="center",zorder=12,bbox=dict(fc="white",ec=COE,lw=0.8,alpha=0.9,pad=0.25))
def arrows(ax,vw,win=None,clen=1.0):
    o=R[vw]["out"]["del"]; X=R[vw]["X"]; v=R[vw]["v"]; F=o["F"]; mv=o["mv"]; sl=np.abs(X[:,2])<1.5
    if win is not None: sl&=(np.abs(X[:,0]-win[0])<win[2]+3)&(np.abs(X[:,1]-win[1])<win[2]+3)
    lf=np.log10(np.maximum(o["Fm"],1e-300)); q=np.zeros(len(X))
    for t in (True,False):
        m=v==t; q[m]=np.argsort(np.argsort(lf[m]))/max(m.sum()-1,1)
    D=np.where(v,2*R[vw]["r"],clen)*(0.6+0.4*q); u=F[:,:2]/np.maximum(np.linalg.norm(F[:,:2],axis=1),1e-300)[:,None]
    for m,c,w,z in ((sl&v&~mv,"#2e86c1",0.0075,3),(sl&~v&~mv,"#c0392b",0.0045,4),(sl&v&mv,EB.HL,0.0075,5),(sl&~v&mv,EB.HL,0.0045,6)):
        ax.quiver(X[m,0],X[m,1],(u[m,0]*D[m]),(u[m,1]*D[m]),color=c,angles="xy",scale_units="xy",scale=1,pivot="middle",width=w*(1.8 if win else 1),
                  headwidth=3.2,headlength=3.6,headaxislength=3.2,zorder=z)
def panel(ax,vw,win=None):
    cp=Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="none",zorder=0); ax.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=ax.transData))
    ax.add_patch(Circle((0,0),RMIN,fill=False,ec="k",lw=2.0,zorder=7)); ax.add_patch(Circle((0,0),30,fill=False,ec="0.4",lw=1.1,ls=":",zorder=7)); ax.plot(0,0,"k+",ms=11,mew=2,zorder=13)
    arrows(ax,vw,win,clen=1.0 if win is None else 0.8)
    if win is None: coe(ax); ax.set_xlim(-32,32); ax.set_ylim(-32,32); ax.set_xticks([-30,0,30]); ax.set_yticks([-30,0,30])
    else: ax.set_xlim(win[0]-win[2],win[0]+win[2]); ax.set_ylim(win[1]-win[2],win[1]+win[2]); ax.set_xlabel("x (Gpc)",fontsize=EB.FS); ax.set_ylabel("y (Gpc)",fontsize=EB.FS)
    ax.set_aspect("equal"); ax.tick_params(labelsize=EB.FK)
def legend(fig,zoom=False):
    hd=[Line2D([],[],color="#c0392b",lw=2,marker=">",ms=6,label="cluster: arrow along its net force"),Line2D([],[],color="#2e86c1",lw=3,marker=">",ms=8,label="void: arrow along its net force (~3 Gpc)"),
        Line2D([],[],color=EB.HL,lw=3,marker=">",ms=8,label="moved body (old spot BEFORE, new spot AFTER)"),Line2D([],[],color="k",lw=2,label="CMB ring / our sphere")]
    if not zoom: hd+=[Line2D([],[],color=COE,lw=3.6,marker=">",ms=9,label="vector from us toward the CoE"),Rectangle((0,0),1,1,fc="#fadbd8",label="CoE core cap"),Line2D([],[],color="0.4",ls=":",lw=1.2,label="30 Gpc edge of model")]
    hd+=[Rectangle((0,0),1,1,fc=EB.BAND,ec="#58d68d",ls="--",lw=1,label="±1× CMB band (pale mint)")]
    fig.legend(handles=hd,loc="upper center",ncol=3,frameon=False,fontsize=EB.FL,bbox_to_anchor=(0.5,0.03),handlelength=2.0,columnspacing=1.2,
               title="arrow length = body size × (0.6–1.0, by rank of log|F|); direction = in-plane net force (Delaunay neighbours)",title_fontsize=EB.FK)
def figure(path,win=None):
    fig=plt.figure(figsize=(13,10.6)); GS=fig.add_gridspec(2,2,height_ratios=[1.0,0.42],hspace=0.36,wspace=0.12 if win is None else 0.22,left=0.06,right=0.98,top=0.9,bottom=0.08)
    for k,vw in enumerate(("before","after")):
        a_=fig.add_subplot(GS[0,k]); panel(a_,vw,win); st=info["del"][vw]["force"]["moved"]
        a_.set_title(f"{vw.upper()}: moved bodies median |F| {st['median_absF_N']:.2e} N, {100*st['frac_toward_ring']:.0f}% toward ring",fontsize=EB.FT)
    EB.residual_below(fig,GS[1,:],d,LAY["before"][4])
    t="foam s2 — net pair force on each body (slab |z| < 1.5 Gpc), BEFORE | AFTER"+("" if win is None else f"\nzoom: {2*win[2]:.0f} × {2*win[2]:.0f} Gpc around x = {win[0]:.1f}, y = {win[1]:.1f} Gpc (busiest moved region on the ring)")
    fig.suptitle(t,fontsize=EB.FT+1,fontweight="bold",color="#8e44ad"); legend(fig,zoom=win is not None); fig.savefig(path,dpi=150,bbox_inches="tight"); plt.close(fig)
figure("net_force_arrows.png")
Xa=R["after"]["X"]; m=R["after"]["out"]["del"]["mv"]&(np.abs(Xa[:,2])<1.5); best=None
for th in np.linspace(0,2*np.pi,72,endpoint=False):
    c=RMIN*np.array([np.cos(th),np.sin(th)]); k=int(np.sum(np.all(np.abs(Xa[m,:2]-c)<5,axis=1)))
    if best is None or k>best[0]: best=(k,c)
figure("net_force_arrows_zoom.png",win=(best[1][0],best[1][1],5.0)); print("zoom",best)

"""Side view of the non-crossing link set: thin slab |z|<H through the centre, links with both ends in the slab.
Kept links coloured by type, links removed by the rule drawn grey dashed. Void cores drawn as their slab cross-sections."""
import os,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
ME=os.path.dirname(os.path.abspath(__file__)); d=np.load(os.path.join(ME,"results","planar_links.npz"))
E,T,kept,why,pos,Rc=d["E"],d["T"],d["kept"],d["why"],d["pos"],d["Rv_core"]; ncl=len(d["MC"])
H=300.0; W=2600.0
ins=(np.abs(pos[:,2])<H)&(np.abs(pos[:,0])<W)&(np.abs(pos[:,1])<W)
fig,axs=plt.subplots(1,2,figsize=(17,9.6))
COL={0:("#e74c3c","cluster–cluster, direct"),1:("#e67e22","cluster–cluster, via one empty node"),2:("#8e44ad","cluster–void"),3:("#16a085","void–void")}
for ax,(lo,hi) in zip(axs,((-W,W),(-1100,1100))):
    for j in np.where(ins[ncl:])[0]+ncl:
        r2=Rc[j]**2-pos[j,2]**2
        if r2>0: ax.add_patch(Circle(pos[j,:2],np.sqrt(r2),fc="#d6eaf8",ec="#5dade2",lw=.6,zorder=0))
    sel=ins[E[:,0]]&ins[E[:,1]]
    for t,(c,lab) in COL.items():
        m=sel&(T==t)&kept; segs=np.stack([pos[E[m,0],:2],pos[E[m,1],:2]],1)
        for s in segs: ax.plot(s[:,0],s[:,1],color=c,lw=1.1,zorder=2)
        ax.plot([],[],color=c,lw=1.6,label=f"kept: {lab} ({m.sum()} in slab)")
    m=sel&~kept
    for s in np.stack([pos[E[m,0],:2],pos[E[m,1],:2]],1): ax.plot(s[:,0],s[:,1],color="0.35",lw=.9,ls=(0,(3,2)),zorder=1)
    ax.plot([],[],color="0.35",ls=(0,(3,2)),label=f"removed by the no-crossing rule ({m.sum()} in slab)")
    c=np.where(ins[:ncl])[0]; ax.scatter(pos[c,0],pos[c,1],s=14,c="k",zorder=3,label="cluster (>3e15 Msun)")
    v=np.where(ins[ncl:])[0]+ncl; ax.scatter(pos[v,0],pos[v,1],s=10,marker="x",c="#2471a3",zorder=3,label="void centre")
    ax.add_patch(Circle((0,0),0,fc="none")); ax.plot([],[],marker="o",ms=10,mfc="#d6eaf8",mec="#5dade2",ls="",label="void core (slab cut)")
    ax.plot(0,0,"*",ms=14,c="gold",mec="k",zorder=4,label="us")
    ax.set_xlim(lo,hi); ax.set_ylim(lo,hi); ax.set_aspect("equal"); ax.set_xlabel("x (Mpc)"); ax.set_ylabel("y (Mpc)")
axs[0].set_title(f"Non-crossing network, slab |z| < {H:.0f} Mpc (links with both ends in the slab)"); axs[1].set_title("Zoom on the centre")
h,l=axs[0].get_legend_handles_labels(); fig.legend(h,l,loc="lower center",ncol=3,fontsize=10,frameon=False)
fig.text(0.5,0.135,"Lines that seem to cross here are a projection effect: in 3D they pass each other more than 2 tube radii (48 Mpc) apart.",ha="center",fontsize=10.5)
fig.tight_layout(rect=(0,0.15,1,1)); fig.savefig(os.path.join(ME,"side_view_noncrossing.png"),dpi=100)

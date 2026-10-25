"""Plain side view of one generated network: slab |z| < half-thickness, void core spheres cut by the slab,
clusters, and every void-cluster link with both ends in the slab. Full region + zoom."""
import os,sys,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection,PatchCollection
ME='/workspace/coe-cosmology-gh/draft/current/analyses/void_mass_network'; sys.path.insert(0,ME); import vmn_gen as G
Mmin=float(sys.argv[1]) if len(sys.argv)>1 else 3e15
g=G.generate(Mmin=Mmin,seed=2000,comp_mode="lcdm"); S=g["S"]/1000; XC=g["XC"]/1000; Rv=g["Rv"]/1000; L=g["links"]; R=g["params"]["R_out"]/1000
def draw(ax,lim,hz,cen=(0,0)):
    sv=(np.abs(S[:,2])<hz)&(np.abs(S[:,0]-cen[0])<lim+2)&(np.abs(S[:,1]-cen[1])<lim+2)
    rc=np.sqrt(np.maximum(Rv[sv]**2-np.maximum(np.abs(S[sv,2])-0,0)**2,0))      # radius of the slice through z=0
    pc=PatchCollection([plt.Circle(S[i,:2],r) for i,r in zip(np.where(sv)[0],rc) if r>0],facecolor="#2471a3",alpha=0.35,edgecolor="none",zorder=1)
    ax.add_collection(pc)
    sc=(np.abs(XC[:,2])<hz)&(np.abs(XC[:,0]-cen[0])<lim+2)&(np.abs(XC[:,1]-cen[1])<lim+2)
    lk=L[sc[L[:,0]]&sv[L[:,1]]]
    seg=np.stack([XC[lk[:,0],:2],S[lk[:,1],:2]],1)
    ax.add_collection(LineCollection(seg,colors="0.3",linewidths=0.35 if lim>8 else 0.8,alpha=0.6,zorder=2))
    m=0.01*g["MC"][sc]; ax.scatter(XC[sc,0],XC[sc,1],s=(2 if lim>8 else 8)+(30 if lim>8 else 80)*(m/m.max())**1.5,c="#c0392b",edgecolors="none",zorder=3); ax.text(0.02,0.02,f"cluster masses at 1%%: {m.min():.2e} to {m.max():.2e} Msun (dot size relative)",transform=ax.transAxes,fontsize=10,bbox=dict(fc="w"))
    ax.add_patch(plt.Circle((0,0),14.3,fill=False,color="k",lw=1.4,zorder=4)); ax.plot(0,0,"k+",ms=12,mew=2,zorder=5)
    ax.add_patch(plt.Circle((0,0),0.5,fill=False,color="k",lw=0.8,ls=":",zorder=4))
    for rr,cc in ((0.23,"m"),(1,"tab:orange"),(2,"g"),(3,"c"),(4,"tab:purple")): ax.add_patch(plt.Circle((0,0),rr,fill=False,color=cc,lw=1.2,ls="--",zorder=4))
    ax.set_xlim(cen[0]-lim,cen[0]+lim); ax.set_ylim(cen[1]-lim,cen[1]+lim); ax.set_aspect("equal"); ax.set_xlabel("Gpc"); ax.set_ylabel("Gpc")
    return sv.sum(),sc.sum(),len(lk)
fig,ax=plt.subplots(2,1,figsize=(10,20.5))
a,b,c=draw(ax[0],R+0.5,0.4); ax[0].set_title(f"Web used on red (seed 2000, lcdm) at 1% mass (full-mass cutoff {Mmin:.0e} Msun, i.e. {0.01*Mmin:.0e} at 1%): slab |z| < 0.4 Gpc of the {R:.0f} Gpc region\n{b} clusters, {c} void–cluster links shown; black circle = our 14.3 Gpc sphere, + = us",fontsize=12)
a,b,c=draw(ax[1],2.5,0.4,cen=(0,0)); ax[1].set_title(f"Zoom (5 × 5 Gpc) around us; dashed = 0.23/1/2/3/4 Gpc evaluation radii: blue = void cores (non-overlapping),\nred = clusters (size ~ mass), grey = links from each cluster to every void it borders",fontsize=12)
fig.text(0.5,0.008,"Blue discs = core spheres of voids centred in the slab (each void is the whole cell between its clusters).\nRed dots = clusters at web nodes. Grey lines = bipartite links (void ↔ cluster only),\ndrawn when both ends are in the slab. Dotted circle = 0.5 Gpc zone kept empty around us.",ha="center",fontsize=10.5)
fig.tight_layout(rect=(0,0.045,1,1)); fig.savefig("/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/constrained/red_web_s2000_bodies_1pct.png",dpi=110)

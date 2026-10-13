"""Full triangulation views for one minimal-edit run (Chris, 10 Oct): Delaunay triangulation of the side-view slab nodes (|z|<1.5 Gpc,
in-plane projection), BEFORE and AFTER separately. Every node is connected; only triangles with an edge longer than LMAX Gpc are dropped
(they span empty gaps / the hull). Figures: tri_<pre>_bodies_mesh.png (bodies + mesh), tri_<pre>_mesh_before_after.png (bare mesh),
tri_<pre>_area_change.png (triangle area, AFTER - BEFORE, on a fixed grid). CoE direction = vector from us (frame30 default)."""
import sys,json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import SymLogNorm
from matplotlib.path import Path
from matplotlib.patches import Circle,Rectangle
from matplotlib.lines import Line2D
from scipy.spatial import Delaunay
import edit_before_after as EB
from common import *
import os
LMAX=float(os.environ.get("TRI_LMAX","inf"))     # full triangulation by default (Chris, 10 Oct): no edge cutoff, no face drops
def tri(X,r):
    sl=np.where((np.abs(X[:,2])<1.5)&(r>0)&(np.linalg.norm(X,axis=1)<=30))[0]; P=X[sl,:2]; S=Delaunay(P).simplices
    if np.isfinite(LMAX): L=np.max([np.linalg.norm(P[S[:,i]]-P[S[:,(i+1)%3]],axis=1) for i in range(3)],axis=0); S=S[L<=LMAX]
    E=np.unique(np.sort(np.r_[S[:,[0,1]],S[:,[1,2]],S[:,[0,2]]],axis=1),axis=0)
    q=P[S]; A=0.5*np.abs((q[:,1,0]-q[:,0,0])*(q[:,2,1]-q[:,0,1])-(q[:,2,0]-q[:,0,0])*(q[:,1,1]-q[:,0,1]))
    return sl,P,S,E,A
def links(ax,P,E,G,lw=(0.4,0.9)):
    mv=G[E[:,0]]|G[E[:,1]]
    ax.add_collection(LineCollection(np.stack([P[E[~mv,0]],P[E[~mv,1]]],1),colors="0.55",lw=lw[0],zorder=2.5))
    ax.add_collection(LineCollection(np.stack([P[E[mv,0]],P[E[mv,1]]],1),colors=EB.HL,lw=lw[1],zorder=3.5)); return int(mv.sum())
def frame(ax):
    cp=Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="none",zorder=0); ax.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=ax.transData))
    ax.add_patch(Circle((0,0),RMIN,fill=False,ec="#1f4e79",lw=2.0,zorder=6)); ax.plot(0,0,"k+",ms=11,mew=2,zorder=9); frame30(ax,lw=1.2,fs=EB.FK)
    ax.set_xticks([-30,0,30]); ax.set_yticks([-30,0,30]); ax.tick_params(labelsize=EB.FK)
COMMON=[Line2D([],[],color="#1f4e79",lw=2,label="CMB ring / our sphere"),Rectangle((0,0),1,1,fc="#fadbd8",label="CoE core cap (inside 30 Gpc)"),
        Line2D([],[],color="0.4",ls=":",lw=1.2,label="30 Gpc edge of model"),Line2D([],[],color="#a93226",lw=2,marker=">",label="vector from us toward the CoE")]
js,pre=sys.argv[1],sys.argv[2]; d=json.load(open(js)); base=d["tag"].replace("_entropy","")
Lb=EB.layout(base); La=EB.layout(base,d); idx=np.array(d["idx"],int); HL=np.zeros(len(Lb[0]),bool); HL[idx[np.array(d["moved"],bool)]]=True
T={k:tri(L[0],L[1]) for k,L in (("before",Lb),("after",La))}; info={}
for k,(sl,P,S,E,A) in T.items():
    G=HL[sl]; info[k]=dict(nodes=len(sl),links=len(E),triangles=len(S),moved_nodes=int(G.sum()),links_touching_moved=int((G[E[:,0]]|G[E[:,1]]).sum()),
                           total_area=float(A.sum()),median_area=float(np.median(A)))
TITLE="foam growth, seed 2 — entropy-rule edits, full triangulation of the slab"
# 1) bodies + mesh
fig=plt.figure(figsize=(13,10.6)); GS=fig.add_gridspec(2,2,height_ratios=[1.0,0.42],hspace=0.36,wspace=0.1,left=0.06,right=0.98,top=0.9,bottom=0.08); ax=[fig.add_subplot(GS[0,0]),fig.add_subplot(GS[0,1])]
EB.residual_below(fig,GS[1,:],d,Lb[4])
fig.suptitle(TITLE+" (bodies + mesh)",fontsize=EB.FT+1,fontweight="bold",color="#8e44ad")
for a_,k,L in zip(ax,("before","after"),(Lb,La)):
    X,r,v,M,_=L; EB.draw_map(a_,X,r,v,M,f"{k.upper()}: {info[k]['links']:,} links, {info[k]['triangles']:,} triangles",EB.FT and "#8e44ad",hl=HL)
    for t_ in a_.texts: t_.set_fontsize(EB.FK)
    sl,P,S,E,A=T[k]; links(a_,P,E,HL[sl],lw=(0.3,0.7))
fig.legend(handles=EB.legend_handles(False)[:4]+[Line2D([],[],color="0.55",lw=1,label="triangulation link"),Line2D([],[],color=EB.HL,lw=1.5,label="link touching a moved body")]+COMMON+[Rectangle((0,0),1,1,fc=EB.BAND,ec="#58d68d",ls="--",lw=1,label="±1× CMB band (pale mint)")],
           loc="upper center",ncol=3,frameon=False,fontsize=EB.FL,bbox_to_anchor=(0.5,0.03),handlelength=1.6,columnspacing=1.2)
fig.savefig(f"tri_{pre}_bodies_mesh.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# 2) bare mesh
fig=plt.figure(figsize=(13,10.6)); GS=fig.add_gridspec(2,2,height_ratios=[1.0,0.42],hspace=0.36,wspace=0.1,left=0.06,right=0.98,top=0.9,bottom=0.08); ax=[fig.add_subplot(GS[0,0]),fig.add_subplot(GS[0,1])]
EB.residual_below(fig,GS[1,:],d,Lb[4])
fig.suptitle(TITLE+" (bare mesh)",fontsize=EB.FT+1,fontweight="bold",color="#8e44ad")
for a_,k,L in zip(ax,("before","after"),(Lb,La)):
    sl,P,S,E,A=T[k]; G=HL[sl]; vs=L[2][sl]; frame(a_); links(a_,P,E,G)
    for m_,c,s in ((~vs&~G,"#c0392b",6),(vs&~G,"#2e86c1",10)): a_.scatter(P[m_,0],P[m_,1],s=s,c=c,lw=0,zorder=4)
    a_.scatter(P[G,0],P[G,1],s=10,c=EB.HL,edgecolor=EB.HLE,lw=0.4,zorder=5)
    a_.set_title(f"{k.upper()}: {info[k]['links']:,} links, {info[k]['nodes']:,} nodes\n{info[k]['links_touching_moved']:,} links touch {info[k]['moved_nodes']} moved bodies",fontsize=EB.FT)
fig.legend(handles=[Line2D([],[],color="0.55",lw=1,label="triangulation link"),Line2D([],[],color=EB.HL,lw=1.5,label="link touching a moved body"),
    Line2D([],[],marker="o",ls="",ms=4,color="#c0392b",label="cluster"),Line2D([],[],marker="o",ls="",ms=5,color="#2e86c1",label="void centre"),
    Line2D([],[],marker="o",ls="",ms=5,mfc=EB.HL,mec=EB.HLE,label="moved body (old spot BEFORE, new spot AFTER)")]+COMMON+[Rectangle((0,0),1,1,fc=EB.BAND,ec="#58d68d",ls="--",lw=1,label="±1× CMB band (pale mint)")],
    loc="upper center",ncol=3,frameon=False,fontsize=EB.FL,bbox_to_anchor=(0.5,0.03),handlelength=1.6,columnspacing=1.2,title=f"nodes and links in the slab |z| < 1.5 Gpc; full Delaunay triangulation (no cutoff)",title_fontsize=EB.FK)
fig.savefig(f"tri_{pre}_mesh_before_after.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# 3) area change on a fixed grid
g=np.linspace(-29.5,29.5,296); GX,GY=np.meshgrid(g,g); pts=np.c_[GX.ravel(),GY.ravel()]
def grid(k):
    sl,P,S,E,A=T[k]; out=np.full(len(pts),np.nan); f=Delaunay(P).find_simplex(pts); D=Delaunay(P).simplices
    key={tuple(sorted(s)):a for s,a in zip(S,A)}; ok=f>=0
    out[ok]=[key.get(tuple(sorted(D[j])),np.nan) for j in f[ok]]; return out.reshape(GX.shape)
Gb,Ga=grid("before"),grid("after"); both=np.isfinite(Gb)&np.isfinite(Ga); D=np.where(both,Ga-Gb,np.nan)
dr=np.abs(np.hypot(GX,GY)-RMIN)
cov=lambda Gx,m: float(np.isfinite(Gx[m&(np.hypot(GX,GY)<=29.5)]).mean())
info["grid"]=dict(covered_0_3_before=cov(Gb,dr<3),covered_0_3_after=cov(Ga,dr<3),covered_6up_before=cov(Gb,dr>=6),covered_6up_after=cov(Ga,dr>=6),frac_grew=float((D[both]>0.01).mean()),frac_shrank=float((D[both]<-0.01).mean()),mean_change=float(np.nanmean(D)),
                  ring_0_3_mean_change=float(np.nanmean(D[both&(dr<3)])),ring_3_6_mean_change=float(np.nanmean(D[both&(dr>=3)&(dr<6)])),ring_6up_mean_change=float(np.nanmean(D[both&(dr>=6)])))
fig,a_=plt.subplots(figsize=(8,8.6)); lim=np.nanpercentile(np.abs(D),98)
im=EB.change_map(a_,GX,GY,D,lim,0.05,hatch=False)
a_.add_patch(Circle((0,0),RMIN,fill=False,ec="k",lw=1.8,zorder=6))
for rr in (RMIN-3,RMIN+3): a_.add_patch(Circle((0,0),rr,fill=False,ec="k",lw=0.9,ls=":",zorder=6))
a_.plot(0,0,"k+",ms=11,mew=2,zorder=9); frame30(a_,lw=1.0,fs=EB.FK); a_.set_xlim(-30.5,30.5); a_.set_ylim(-30.5,30.5); a_.set_xticks([-30,0,30]); a_.set_yticks([-30,0,30]); a_.tick_params(labelsize=EB.FK)
a_.set_title(f"Triangle area change, AFTER − BEFORE (foam s2)\n{100*info['grid']['frac_grew']:.0f}% of the area grew, {100*info['grid']['frac_shrank']:.0f}% shrank",fontsize=EB.FT)
cb=fig.colorbar(im,ax=a_,orientation="horizontal",fraction=0.045,pad=0.08,ticks=[-1,-0.1,0,0.1,1]); cb.ax.set_xticklabels(["−1","−0.1","0","0.1","1"]); cb.set_label("Δ triangle area at that point (Gpc²): blue = compressed (smaller), red = expanded (bigger), white = no change",fontsize=EB.FK); cb.ax.tick_params(labelsize=EB.FK-1)
fig.legend(handles=[Line2D([],[],color="k",lw=1.8,label="CMB ring"),Line2D([],[],color="k",ls=":",lw=1,label="±3 Gpc from the ring"),
    Line2D([],[],color="#a93226",lw=2,marker=">",label="vector from us toward the CoE"),Line2D([],[],color="0.4",ls=":",lw=1,label="30 Gpc edge; map clipped to the area inside both views' triangulations")],
    loc="upper center",bbox_to_anchor=(0.5,0.0),ncol=2,frameon=False,fontsize=EB.FL)
fig.savefig(f"tri_{pre}_area_change.png",dpi=140,bbox_inches="tight"); plt.close(fig)
json.dump(info,open(f"tri_{pre}_info.json","w"),indent=1); print(json.dumps(info,indent=1))

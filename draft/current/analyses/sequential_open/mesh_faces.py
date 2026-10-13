"""Face areas of the planar non-crossing mesh (mesh_view.py) for one run, BEFORE vs AFTER (Chris, 10 Oct).
Faces: half-edge walk of the straight-line planar graph (in-plane projection, slab |z|<1.5 Gpc). Each component's outer face (clockwise,
negative signed area) is dropped; faces with any vertex beyond 28 Gpc (cut by the 30 Gpc edge) are dropped; faces larger than AMAX Gpc^2
are counted as slab holes (open gaps where the slab has no nodes) and dropped."""
import sys,json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm,TwoSlopeNorm,SymLogNorm
from matplotlib.path import Path
from matplotlib.lines import Line2D
from matplotlib.patches import Circle,Rectangle
import mesh_view as MV, edit_before_after as EB
from common import *
AMAX=40.0
def faces(P,E):
    n=len(P); nb=[[] for _ in range(n)]
    for a,b in E: nb[a].append(b); nb[b].append(a)
    order={}
    for u in range(n):
        if not nb[u]: continue
        ang=np.arctan2(P[nb[u],1]-P[u,1],P[nb[u],0]-P[u,0]); s=[nb[u][k] for k in np.argsort(ang)]; order[u]={v:i for i,v in enumerate(s)}; nb[u]=s
    seen=set(); F=[]
    for a,b in E:
        for u,v in ((a,b),(b,a)):
            if (u,v) in seen: continue
            f=[]; x,y=u,v
            while (x,y) not in seen:
                seen.add((x,y)); f.append(x); i=order[y][x]; z=nb[y][(i-1)%len(nb[y])]; x,y=y,z     # next edge clockwise at y -> faces CCW
            F.append(f)
    out=[]
    for f in F:
        q=P[f]; A=0.5*np.sum(q[:,0]*np.roll(q[:,1],-1)-np.roll(q[:,0],-1)*q[:,1]); out.append((f,A))
    return out
def valid(P,FA):
    keep=[]; stats=dict(outer=0,edge=0,hole=0)
    for f,A in FA:
        if A<=0: stats["outer"]+=1; continue
        if np.max(np.hypot(P[f,0],P[f,1]))>28: stats["edge"]+=1; continue
        if A>AMAX: stats["hole"]+=1; continue
        keep.append((f,A))
    return keep,stats
js=sys.argv[1]; pre=sys.argv[2]; d=json.load(open(js)); base=d["tag"].replace("_entropy","")
Xb,rb,vb,_,_=EB.layout(base); Xa,ra,va,_,_=EB.layout(base,d); idx=np.array(d["idx"],int); HL=np.zeros(len(Xb),bool); HL[idx[np.array(d["moved"],bool)]]=True
V={}
for k,(X,r,v) in (("before",(Xb,rb,vb)),("after",(Xa,ra,va))):
    sl,E,T=MV.network(X,r,v); P=X[sl,:2]; FA,st=valid(P,faces(P,E)); A=np.array([a for _,a in FA])
    cen=np.array([P[f].mean(0) for f,_ in FA]); dr=np.abs(np.hypot(cen[:,0],cen[:,1])-RMIN)
    V[k]=dict(sl=sl,E=E,P=P,F=FA,A=A,dr=dr,st=st)
# (1) per moved body: summed area of valid faces touching it (only bodies in the slab in both views)
per=[]
for i in np.where(HL)[0]:
    s=[]
    for k in ("before","after"):
        loc=np.where(V[k]["sl"]==i)[0]
        if not len(loc): s=None; break
        j=loc[0]; s.append(sum(a for f,a in V[k]["F"] if j in f))
    if s: per.append((i,s[0],s[1],abs(np.hypot(Xb[i,0],Xb[i,1])-RMIN)))
per=np.array(per)
# (2) grid map
g=np.linspace(-28,28,281); GX,GY=np.meshgrid(g,g); pts=np.c_[GX.ravel(),GY.ravel()]
def grid(k):
    out=np.full(len(pts),np.nan); P=V[k]["P"]
    for f,A in V[k]["F"]:
        q=P[f]; lo,hi=q.min(0),q.max(0); m=np.where(np.all((pts>=lo)&(pts<=hi),axis=1))[0]
        if len(m): out[m[Path(q).contains_points(pts[m])]]=A
    return out.reshape(GX.shape)
Gb,Ga=grid("before"),grid("after"); D=Ga-Gb
# numbers
def band(k,lo,hi): m=(V[k]["dr"]>=lo)&(V[k]["dr"]<hi); return float(np.median(V[k]["A"][m])) if m.any() else None,int(m.sum())
R={k:dict(n_faces=len(V[k]["A"]),total=float(V[k]["A"].sum()),median=float(np.median(V[k]["A"])),mean=float(V[k]["A"].mean()),dropped=V[k]["st"],
          ring_0_3=band(k,0,3),ring_3_6=band(k,3,6),ring_6up=band(k,6,99)) for k in V}
R["per_moved"]=dict(n=len(per),sum_before=float(per[:,1].sum()),sum_after=float(per[:,2].sum()),median_ratio=float(np.median(per[:,2]/np.maximum(per[:,1],1e-9))),
    n_grew=int((per[:,2]>per[:,1]).sum()),n_shrank=int((per[:,2]<per[:,1]).sum()),
    near_ring_0_3=dict(n=int((per[:,3]<3).sum()),sum_before=float(per[per[:,3]<3,1].sum()),sum_after=float(per[per[:,3]<3,2].sum())))
both=np.isfinite(Gb)&np.isfinite(Ga); R["grid"]=dict(points_both=int(both.sum()),frac_grew=float((D[both]>0.01).mean()),frac_shrank=float((D[both]<-0.01).mean()),mean_change=float(D[both].mean()))
dGr=np.abs(np.hypot(GX,GY)-RMIN); m=both&(dGr<3); R["grid"]["ring_0_3_mean_change"]=float(D[m].mean()); R["grid"]["ring_0_3_mean_before"]=float(Gb[m].mean())
json.dump(R,open(pre+"_faces.json","w"),indent=1); print(json.dumps(R,indent=1))
# (a) faces filled by area
norm=LogNorm(vmin=0.05,vmax=AMAX); cm=plt.get_cmap("viridis")
fig,ax=plt.subplots(1,2,figsize=(13,7.4)); fig.subplots_adjust(wspace=0.1,left=0.05,right=0.98,top=0.86,bottom=0.16)
fig.suptitle("foam growth, seed 2 — mesh faces coloured by area (same scale in both views)",fontsize=EB.FT+1,fontweight="bold",color="#8e44ad")
for a_,k in zip(ax,("before","after")):
    P=V[k]["P"]; pc=PolyCollection([P[f] for f,_ in V[k]["F"]],array=V[k]["A"],cmap=cm,norm=norm,edgecolors="0.35",linewidths=0.25,zorder=1); a_.add_collection(pc)
    a_.add_patch(Circle((0,0),RMIN,fill=False,ec="#e74c3c",lw=2.0,zorder=6)); a_.plot(0,0,"k+",ms=11,mew=2,zorder=7); frame30(a_,lw=1.2)
    for t in a_.texts: t.set_fontsize(EB.FK)
    a_.set_xticks([-30,0,30]); a_.set_yticks([-30,0,30]); a_.tick_params(labelsize=EB.FK)
    a_.set_title(f"{k.upper()}: {R[k]['n_faces']} faces, median {R[k]['median']:.2f} Gpc²",fontsize=EB.FT)
cb=fig.colorbar(pc,ax=ax,orientation="horizontal",fraction=0.045,pad=0.06,aspect=50); cb.set_label("face area (Gpc², log scale)",fontsize=EB.FS); cb.ax.tick_params(labelsize=EB.FK)
fig.legend(handles=[Line2D([],[],color="#e74c3c",lw=2,label="CMB ring / our sphere"),Line2D([],[],color="0.4",ls=":",lw=1.2,label="30 Gpc edge of model"),
    Line2D([],[],color="#a93226",marker=">",label="direction to the CoE"),Rectangle((0,0),1,1,fc="white",ec="0.6",label="white = no face (outer face, slab hole, or cut by the edge)")],
    loc="upper center",bbox_to_anchor=(0.5,0.035),ncol=2,frameon=False,fontsize=EB.FL)
fig.savefig(pre+"_faces_area.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# (b) change map
fig,a_=plt.subplots(figsize=(8,8.4)); lim=np.nanpercentile(np.abs(D),98)
im=EB.change_map(a_,GX,GY,np.where(both,D,np.nan),lim,0.1)
a_.add_patch(Circle((0,0),RMIN,fill=False,ec="k",lw=1.8,zorder=6))
for rr in (RMIN-3,RMIN+3): a_.add_patch(Circle((0,0),rr,fill=False,ec="k",lw=0.9,ls=":",zorder=6))
a_.plot(0,0,"k+",ms=11,mew=2); a_.set_xlim(-30,30); a_.set_ylim(-30,30); a_.set_aspect("equal"); a_.set_xticks([-30,0,30]); a_.set_yticks([-30,0,30]); a_.tick_params(labelsize=EB.FK)
frame30(a_,lw=1.0,fs=EB.FK); a_.set_xlim(-30.5,30.5); a_.set_ylim(-30.5,30.5)
a_.set_title(f"Change in local face area, AFTER − BEFORE\n{100*R['grid']['frac_grew']:.0f}% of the area grew, {100*R['grid']['frac_shrank']:.0f}% shrank",fontsize=EB.FT)
cb=fig.colorbar(im,ax=a_,orientation="horizontal",fraction=0.045,pad=0.08); cb.set_label("Δ face area at that point (Gpc²): blue = compressed (smaller), red = expanded (bigger), white = no change",fontsize=EB.FK); cb.ax.tick_params(labelsize=EB.FK-1)
fig.legend(handles=[Line2D([],[],color="k",lw=1.8,label="CMB ring"),Line2D([],[],color="k",ls=":",lw=1,label="±3 Gpc from the ring"),Line2D([],[],color="#a93226",lw=2,marker=">",label="vector from us toward the CoE"),Rectangle((0,0),1,1,**EB.NOFACE,label="light grey = no face in one of the views")],
    loc="upper center",bbox_to_anchor=(0.5,0.0),ncol=2,frameon=False,fontsize=EB.FL)
fig.savefig(pre+"_faces_change.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# (c) distributions
fig,ax=plt.subplots(1,2,figsize=(13,5.6)); fig.subplots_adjust(wspace=0.28,bottom=0.25)
for k,c,ls in (("before","0.45","-"),("after",EB.HL,"--")):
    A=np.sort(V[k]["A"]); ax[0].plot(A,np.arange(1,len(A)+1)/len(A),color=c,ls=ls,lw=2.2,label=f"{k}: {len(A)} faces, total {A.sum():.0f} Gpc², median {np.median(A):.2f} Gpc²")
ax[0].set_xscale("log"); ax[0].set_xlabel("face area (Gpc²)",fontsize=EB.FS); ax[0].set_ylabel("fraction of faces ≤ area",fontsize=EB.FS); ax[0].set_title("All faces (CDF)",fontsize=EB.FT); ax[0].grid(alpha=.3)
lab=["0–3 Gpc","3–6 Gpc","6+ Gpc"]; x=np.arange(3); keys=("ring_0_3","ring_3_6","ring_6up")
for o,k,c in ((-0.21,"before","0.6"),(0.21,"after",EB.HL)):
    vals=[R[k][kk][0] for kk in keys]; ax[1].bar(x+o,vals,0.36,color=c,label=k)
    for xi,vv,kk in zip(x,vals,keys): ax[1].text(xi+o,vv,f"{vv:.2f}\nn={R[k][kk][1]}",ha="center",va="bottom",fontsize=EB.FK-2)
ax[1].set_xticks(x); ax[1].set_xticklabels(lab,fontsize=EB.FK); ax[1].set_xlabel("distance of face centre from the CMB ring",fontsize=EB.FS); ax[1].set_ylabel("median face area (Gpc²)",fontsize=EB.FS)
ax[1].set_title("Median face area by distance from the ring",fontsize=EB.FT); ax[1].set_ylim(0,max(R[k][kk][0] for k in ("before","after") for kk in keys)*1.35)
for a_ in ax: a_.tick_params(labelsize=EB.FK)
h0,l0=ax[0].get_legend_handles_labels(); fig.legend(h0,l0,loc="upper center",bbox_to_anchor=(0.5,0.1),ncol=1,frameon=False,fontsize=EB.FL)
fig.savefig(pre+"_faces_dist.png",dpi=140,bbox_inches="tight")

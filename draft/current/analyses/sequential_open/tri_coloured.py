"""tri_foam_s2_mesh_coloured.png: full Delaunay mesh BEFORE (plain) | AFTER (each triangle filled by A_after - A_before(local), where
A_before(local) = area of the BEFORE triangle containing the AFTER triangle's centroid). Blue = compressed, red = expanded, white = no change."""
import json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.colors import SymLogNorm
from scipy.spatial import Delaunay
src=open("tri_view.py").read(); exec(src[:src.index("js,pre=sys.argv")])
d=json.load(open("minedit_ent/b_foam_rim_s2_tidalcore_entropy.json")); base=d["tag"].replace("_entropy","")
Lb=EB.layout(base); La=EB.layout(base,d); idx=np.array(d["idx"],int); HL=np.zeros(len(Lb[0]),bool); HL[idx[np.array(d["moved"],bool)]]=True
slb,Pb,Sb,Eb,Ab=tri(Lb[0],Lb[1]); sla,Pa,Sa,Ea,Aa=tri(La[0],La[1])
cen=Pa[Sa].mean(1); DT=Delaunay(Pb); j=DT.find_simplex(cen); key={tuple(sorted(s)):a for s,a in zip(Sb,Ab)}
Aloc=np.array([key[tuple(sorted(DT.simplices[k]))] if k>=0 else np.nan for k in j]); dA=Aa-Aloc
lim=np.nanpercentile(np.abs(dA),98); norm=SymLogNorm(linthresh=0.05,vmin=-lim,vmax=lim)
fig=plt.figure(figsize=(13,11.4)); GS=fig.add_gridspec(3,2,height_ratios=[1.0,0.035,0.42],hspace=0.42,wspace=0.1,left=0.06,right=0.98,top=0.9,bottom=0.07)
a0=fig.add_subplot(GS[0,0]); a1=fig.add_subplot(GS[0,1])
for a_,P,E,sl,L,t in ((a0,Pb,Eb,slb,Lb,f"BEFORE: {len(Sb):,} triangles, median {np.median(Ab):.2f} Gpc²"),(a1,Pa,Ea,sla,La,f"AFTER: {len(Sa):,} triangles, median {np.median(Aa):.2f} Gpc²\nfilled by area change vs local BEFORE area")):
    cp=Circle(C_G[:2],A_IN_G,fc="#f1948a",ec="none",alpha=0.25,zorder=1.5); a_.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=a_.transData))
    if a_ is a1:
        pc=PolyCollection(Pa[Sa],array=dA,cmap="RdBu_r",norm=norm,edgecolors="none",zorder=1); a_.add_collection(pc)
    a_.add_collection(LineCollection(np.stack([P[E[:,0]],P[E[:,1]]],1),colors="0.45",lw=0.3,zorder=2))
    vs=L[2][sl]; a_.scatter(P[~vs,0],P[~vs,1],s=3,c="#c0392b",lw=0,zorder=3); a_.scatter(P[vs,0],P[vs,1],s=5,c="#2e86c1",lw=0,zorder=3)
    a_.add_patch(Circle((0,0),RMIN,fill=False,ec="k",lw=2.0,zorder=6)); a_.plot(0,0,"k+",ms=11,mew=2,zorder=9); frame30(a_,lw=1.2,fs=EB.FK)
    a_.set_xticks([-30,0,30]); a_.set_yticks([-30,0,30]); a_.tick_params(labelsize=EB.FK); a_.set_title(t,fontsize=EB.FT)
cax=fig.add_subplot(GS[1,:]); cb=fig.colorbar(pc,cax=cax,orientation="horizontal",ticks=[-1,-0.1,0,0.1,1]); cb.ax.set_xticklabels(["−1","−0.1","0","0.1","1"])
cb.set_label("AFTER triangle area − local BEFORE area (Gpc²): "+EB.CHANGE_LABEL,fontsize=EB.FK); cb.ax.tick_params(labelsize=EB.FK)
EB.residual_below(fig,GS[2,:],d,Lb[4])
fig.suptitle("foam growth, seed 2 — full Delaunay triangulation, before | after coloured by area change",fontsize=EB.FT+1,fontweight="bold",color="#8e44ad")
fig.legend(handles=[Line2D([],[],color="0.45",lw=1,label="triangulation link"),Line2D([],[],marker="o",ls="",ms=4,color="#c0392b",label="cluster"),
    Line2D([],[],marker="o",ls="",ms=5,color="#2e86c1",label="void centre"),Line2D([],[],color="k",lw=2,label="CMB ring / our sphere"),
    Rectangle((0,0),1,1,fc="#fadbd8",label="CoE core cap (inside 30 Gpc)"),Line2D([],[],color="#a93226",lw=2,marker=">",label="vector from us toward the CoE"),
    Rectangle((0,0),1,1,fc=EB.BAND,ec="#58d68d",ls="--",lw=1,label="±1× CMB band (pale mint)")],loc="upper center",ncol=4,frameon=False,fontsize=EB.FL,bbox_to_anchor=(0.5,0.03))
fig.savefig("tri_foam_s2_mesh_coloured.png",dpi=140,bbox_inches="tight")
dr=np.abs(np.hypot(cen[:,0],cen[:,1])-RMIN)
for lab,m in (("whole",np.ones(len(dr),bool)),("0-3",dr<3),("3-6",(dr>=3)&(dr<6)),("6+",dr>=6)):
    mm=m&np.isfinite(dA); print(lab,"n",mm.sum(),"bigger",int((dA[mm]>0.01).sum()),"smaller",int((dA[mm]<-0.01).sum()),"median dA %.3f"%np.median(dA[mm]))

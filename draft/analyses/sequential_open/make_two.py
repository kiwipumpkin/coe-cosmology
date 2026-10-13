"""Priority deliverable: B random (s2) and B foam rule (s3), open region vs old runs."""
import os,json,glob,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from common import *
plt.rcParams.update({"font.size":8,"savefig.dpi":300,"pdf.fonttype":42})
RUN={"seqB_random":("random",2,"#7f8c8d","B: growth, random"),"seqB_foam_rim":("foam_rim",3,"#5dade2","B: growth, foam rule")}
D={}
for q,(rule,seed,col,nm) in RUN.items():
    old=json.load(open(f"{SEQ}/b_{rule}_Rv1.5_s{seed}.json")); new=json.load(open(f"{OUT}/b_{rule}_s{seed}.json")); f=np.load(f"{OUT}/final_b_{rule}_s{seed}.npz")
    D[q]=dict(old=old,new=new,f=f)
fig,ax=plt.subplots(figsize=(3.6,5.4))
for q,(rule,seed,col,nm) in RUN.items():
    o,n=D[q]["old"],D[q]["new"]
    ax.plot(o["traj_obj"],o["traj_left"],color=col,lw=1.0,ls=(0,(4,2)))
    ax.plot(n["traj_obj"],n["traj_left"],color=col,lw=1.5)
ax.axhline(1,color="green",lw=1.5); ax.axhline(237.3,color="0.5",ls=":",lw=0.8); ax.axhline(88.8,color="0.6",ls=":",lw=0.8)
ax.text(1.15,1.12,"CMB flat, 1×",color="green",fontsize=6.5,va="bottom"); ax.text(1.15,250,"start 237×",color="0.4",fontsize=6.5,va="bottom"); ax.text(1.15,93,"four clumps 89×",color="0.45",fontsize=6.5,va="bottom")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1,2e5); ax.set_ylim(0.5,3000)
ax.set_xlabel("bodies placed"); ax.set_ylabel(r"leftover ($\times$ floor)"); ax.set_title("B random and B foam rule:\ncore and our sphere open",fontsize=8)
hd=[Line2D([],[],color=RUN[q][2],lw=1.5,label=RUN[q][3]) for q in RUN]+[Line2D([],[],color="0.3",lw=1.5,label="open region (this run)"),Line2D([],[],color="0.3",lw=1.0,ls=(0,(4,2)),label="old run (14.3–42.9 Gpc, core blocked)")]
ax.legend(handles=hd,loc="upper center",bbox_to_anchor=(0.5,-0.13),frameon=False,fontsize=6.5)
fig.savefig(f"{OUT}/fig7_open_random_foam.png",bbox_inches="tight",dpi=200); fig.savefig(f"{OUT}/fig7_open_random_foam.pdf",bbox_inches="tight"); plt.close(fig)
order=sorted(RUN,key=lambda q:-D[q]["new"]["final_left"])
fig=plt.figure(figsize=(7.1,4.6)); GO=fig.add_gridspec(1,2,wspace=0.25,left=0.07,right=0.98,top=0.92,bottom=0.2)
for i,q in enumerate(order):
    rule,seed,col,nm=RUN[q]; f=D[q]["f"]; v=f["kp"]=="void"; Xv,rv,Xc=f["Xp"][v],f["rp"][v],f["Xp"][~v]
    sub=GO[0,i].subgridspec(2,1,height_ratios=[1.7,0.95],hspace=0.6); ax=fig.add_subplot(sub[0])
    ax.add_patch(Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="#e6b0aa",lw=0.6,ls="--",zorder=0)); ax.text(C_G[0],C_G[1]-A_IN_G+4,"CoE core (open)",color="#a93226",fontsize=6,ha="center",zorder=6,bbox=dict(fc="white",ec="none",alpha=0.8,pad=0.5))
    ax.add_patch(Circle((0,0),42.9,fill=False,ec="0.6",lw=0.6,ls=":"))
    sv=np.abs(Xv[:,2])<1.5; sc=np.abs(Xc[:,2])<1.5
    ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(Xv[sv],rv[sv])],fc="#aed6f1",ec="#2471a3",lw=0.25,zorder=2))
    ax.scatter(Xc[sc,0],Xc[sc,1],s=1.2,c="#c0392b",lw=0,zorder=3)
    ax.add_patch(Circle((0,0),RMIN,fill=False,ec="#1f4e79",lw=1.4,zorder=4)); ax.plot(0,0,"k+",ms=5,mew=1,zorder=5)
    ax.set_xlim(-106,48); ax.set_ylim(-64,48); ax.set_aspect("equal"); ax.tick_params(labelsize=6); ax.set_xlabel("x (Gpc)",fontsize=6.5); ax.set_ylabel("y (Gpc)",fontsize=6.5)
    for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(2.2)
    ax.set_title(f"{i+1}. {nm}\n{v.sum():,} voids + {(~v).sum():,} clusters (slab |z|<1.5 Gpc)",color=col,fontsize=7.2,fontweight="bold")
    r_=fig.add_subplot(sub[1]); res=(R-f["w"])/FLOOR; r_.axhspan(-1,1,color="green",alpha=0.25,lw=0); r_.axhline(1,color="green",lw=0.7); r_.axhline(-1,color="green",lw=0.7); r_.axhline(0,color="0.6",lw=0.5); r_.plot(np.arange(360),res,color=col,lw=1.1)
    lim=1.15*max(np.abs(res).max(),1); r_.set_ylim(-lim,lim); r_.set_xlim(0,360); r_.set_xticks([0,90,180,270,360]); r_.tick_params(labelsize=6)
    r_.set_title(f"leftover {float(rms(R-f['w'])/FLOOR):,.1f}× the CMB level",fontsize=7.2,color=col,fontweight="bold"); r_.set_xlabel("longitude ℓ (°) around the ring",fontsize=6.3); r_.set_ylabel("× $10^{-5}$",fontsize=6.5)
hd=[Line2D([],[],marker="o",ls="",color="#c0392b",ms=3.5,label="cluster"),Line2D([],[],marker="o",ls="",mfc="#aed6f1",mec="#2471a3",ms=7,label="void (to scale)"),
    Line2D([],[],color="#1f4e79",lw=1.4,label="CMB ring / our sphere (open inside)"),plt.Rectangle((0,0),1,1,fc="#fadbd8",ec="#e6b0aa",ls="--",label="CoE core (open to bodies)"),
    Line2D([],[],color="0.6",lw=0.6,ls=":",label="reach outside the core (42.9 Gpc)"),plt.Rectangle((0,0),1,1,fc="green",alpha=0.25,ec="none",label="±1× (CMB level)")]
fig.legend(handles=hd,loc="lower center",ncol=3,frameon=False,fontsize=6.6,bbox_to_anchor=(0.5,0.0))
fig.savefig(f"{OUT}/fig8_open_random_foam.png",bbox_inches="tight",dpi=200); fig.savefig(f"{OUT}/fig8_open_random_foam.pdf",bbox_inches="tight"); plt.close(fig)
for q in RUN:
    o,n=D[q]["old"],D[q]["new"]; print(q,"old final %.0f min %.1f n %d | new final %.0f min %.1f at %d n %d vf %.3f gap %.3f"%(o["final_left"],o["min_left"],o["n_obj"],n["final_left"],n["min_left"],n["argmin"],n["n_obj"],n["void_vol_frac"],n["min_gap_gpc"]),n["where_voids"],n["where_clusters"],"amps",np.round(n["final_amps"],0))

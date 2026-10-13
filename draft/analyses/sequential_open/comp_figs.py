"""Figures for the compensated-web run (Chris, 10 Oct): usage python3 comp_figs.py minedit_comp/comp_sN_entropy.json OUTDIR
comp_web_side_before_after.png, comp_web_mesh_before_after.png, comp_web_hist.png (moderate fonts, legends below, residual below the views)."""
import sys,json,os,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.colors import SymLogNorm
from scipy.spatial import Delaunay
JS,OD=sys.argv[1],sys.argv[2]; os.makedirs(OD,exist_ok=True)
src=open("tri_view.py").read(); exec(src[:src.index("js,pre=sys.argv")])
d=json.load(open(JS)); base=d["tag"].replace("_entropy",""); seed=base.split("_s")[-1]
Lb=EB.layout(base); La=EB.layout(base,d); idx=np.array(d["idx"],int); HL=np.zeros(len(Lb[0]),bool); HL[idx[np.array(d["moved"],bool)]]=True
TITLE=f"Compensated random web, seed {seed} — move + resize, entropy rule"
f2=lambda s_: f"m1–12 {np.hypot(s_['rms_m1_5'],s_['rms_m6_12']):.2f}×, m13+ {s_['rms_13up']:.2f}×, symmetric {s_['sym']*100:.0f}%"
# 1) side view, bodies to scale
fig=plt.figure(figsize=(12,11.2)); G=fig.add_gridspec(2,2,height_ratios=[1.0,0.45],hspace=0.42,wspace=0.12,left=0.08,right=0.98,top=0.88,bottom=0.08)
fig.suptitle(TITLE,fontsize=EB.FT+2,fontweight="bold",y=0.975)
sb=d["before"]; ed=f"{d['n_moved']} moved, {d['n_resized']} resized"
EB.draw_map(fig.add_subplot(G[0,0]),*Lb[:4],f"BEFORE: random start, {len(Lb[0]):,} bodies\nleftover {sb['total']:.0f}×","0.3",hl=HL)
EB.draw_map(fig.add_subplot(G[0,1]),*La[:4],f"AFTER: {d['n_edited']} edited ({100*d['frac_edited']:.1f}%), {ed}\nleftover {d['left']:.2f}×","0.3",hl=HL)
EB.residual_below(fig,G[1,:],d,Lb[4])
fig.legend(handles=EB.legend_handles(False),loc="upper center",ncol=3,frameon=False,fontsize=EB.FL,bbox_to_anchor=(0.5,0.0),handlelength=1.6,columnspacing=1.2,
           title="clusters red, voids blue (drawn to size, slab |z| < 1.5 Gpc); moved bodies green in both views (all z); resized-only keep red/blue",title_fontsize=EB.FK)
fig.savefig(f"{OD}/comp_web_side_before_after.png",dpi=130,bbox_inches="tight"); plt.close(fig)
# 2) mesh: BEFORE plain, AFTER triangles coloured by area change vs the local BEFORE triangle
slb,Pb,Sb,Eb,Ab=tri(Lb[0],Lb[1]); sla,Pa,Sa,Ea,Aa=tri(La[0],La[1])
cen=Pa[Sa].mean(1); DT=Delaunay(Pb); j=DT.find_simplex(cen); key={tuple(sorted(s)):a for s,a in zip(Sb,Ab)}
Aloc=np.array([key[tuple(sorted(DT.simplices[k]))] if k>=0 else np.nan for k in j]); dA=Aa-Aloc
lim=np.nanpercentile(np.abs(dA),98); norm=SymLogNorm(linthresh=0.05,vmin=-lim,vmax=lim)
fig=plt.figure(figsize=(13,11.6)); GS=fig.add_gridspec(3,2,height_ratios=[1.0,0.035,0.42],hspace=0.42,wspace=0.1,left=0.06,right=0.98,top=0.87,bottom=0.07)
a0=fig.add_subplot(GS[0,0]); a1=fig.add_subplot(GS[0,1])
for a_,P,E,sl,L,t in ((a0,Pb,Eb,slb,Lb,f"BEFORE: {len(Sb):,} triangles, median {np.median(Ab):.2f} Gpc²"),(a1,Pa,Ea,sla,La,f"AFTER: {len(Sa):,} triangles, median {np.median(Aa):.2f} Gpc²\nfilled by area change vs local BEFORE area")):
    cp=Circle(C_G[:2],A_IN_G,fc="#f1948a",ec="none",alpha=0.25,zorder=1.5); a_.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=a_.transData))
    if a_ is a1: pc=PolyCollection(Pa[Sa],array=dA,cmap="RdBu_r",norm=norm,edgecolors="none",zorder=1); a_.add_collection(pc)
    a_.add_collection(LineCollection(np.stack([P[E[:,0]],P[E[:,1]]],1),colors="0.45",lw=0.3,zorder=2))
    vs=L[2][sl]; G_=HL[sl]
    a_.scatter(P[~vs&~G_,0],P[~vs&~G_,1],s=3,c="#c0392b",lw=0,zorder=3); a_.scatter(P[vs&~G_,0],P[vs&~G_,1],s=6,c="#2e86c1",lw=0,zorder=3)
    a_.scatter(P[G_,0],P[G_,1],s=8,c=EB.HL,edgecolor=EB.HLE,lw=0.3,zorder=4)
    a_.add_patch(Circle((0,0),RMIN,fill=False,ec="k",lw=2.0,zorder=6)); a_.plot(0,0,"k+",ms=11,mew=2,zorder=9); frame30(a_,lw=1.2,fs=EB.FK)
    a_.set_xticks([-30,0,30]); a_.set_yticks([-30,0,30]); a_.tick_params(labelsize=EB.FK); a_.set_title(t,fontsize=EB.FT)
cax=fig.add_subplot(GS[1,:]); cb=fig.colorbar(pc,cax=cax,orientation="horizontal"); cb.set_label("AFTER triangle area − local BEFORE area (Gpc²): "+EB.CHANGE_LABEL,fontsize=EB.FK); cb.ax.tick_params(labelsize=EB.FK)
EB.residual_below(fig,GS[2,:],d,Lb[4])
fig.suptitle(TITLE+"\nfull Delaunay triangulation of the slab |z| < 1.5 Gpc (no cutoff)",fontsize=EB.FT+1,fontweight="bold")
fig.legend(handles=[Line2D([],[],color="0.45",lw=1,label="triangulation link"),Line2D([],[],marker="o",ls="",ms=4,color="#c0392b",label="cluster"),
    Line2D([],[],marker="o",ls="",ms=5,color="#2e86c1",label="void centre"),Line2D([],[],marker="o",ls="",ms=5,mfc=EB.HL,mec=EB.HLE,label="moved body"),
    Line2D([],[],color="k",lw=2,label="CMB ring / our sphere"),Rectangle((0,0),1,1,fc="#fadbd8",label="CoE core cap (inside 30 Gpc)"),
    Line2D([],[],color="#a93226",lw=2,marker=">",label="vector from us toward the CoE"),Rectangle((0,0),1,1,fc=EB.BAND,ec="#58d68d",ls="--",lw=1,label="±1× CMB band (pale mint)")],
    loc="upper center",ncol=4,frameon=False,fontsize=EB.FL,bbox_to_anchor=(0.5,0.03))
fig.savefig(f"{OD}/comp_web_mesh_before_after.png",dpi=130,bbox_inches="tight"); plt.close(fig)
# 3) histograms, stacked
def stats(L):
    sl,P,S,E,A=tri(L[0],L[1]); c=P[S].mean(1); dr=np.abs(np.hypot(c[:,0],c[:,1])-RMIN)
    deg=np.bincount(E.ravel(),minlength=len(P)); Ll=np.linalg.norm(P[E[:,0]]-P[E[:,1]],axis=1); return A,dr,deg,Ll
B=stats(Lb); Aft=stats(La)
fig,ax=plt.subplots(4,1,figsize=(9,16)); fig.subplots_adjust(hspace=0.5,top=0.92,bottom=0.08)
ab=np.logspace(-3,np.log10(max(B[0].max(),Aft[0].max())*1.05),40); db=np.arange(-0.5,max(B[2].max(),Aft[2].max())+1.5,1); lb=np.linspace(0,max(B[3].max(),Aft[3].max())*1.02,50)
info={}
for i,(t,xb,get,xl,lg) in enumerate((("(a) triangle area, whole slab",ab,lambda S: S[0],"triangle area (Gpc², log bins)",True),
                                       ("(b) triangle area, 0–3 Gpc from the CMB ring",ab,lambda S: S[0][S[1]<3],"triangle area (Gpc², log bins)",True),
                                       ("(c) connections per node",db,lambda S: S[2],"connections per node",False),("(d) link length",lb,lambda S: S[3],"link length (Gpc)",False))):
    a_=ax[i]; x0=get(B); x1=get(Aft); info[t]=dict(before_n=int(len(x0)),after_n=int(len(x1)),before_median=float(np.median(x0)),after_median=float(np.median(x1)),before_mean=float(np.mean(x0)),after_mean=float(np.mean(x1)))
    a_.hist(x0,xb,color="0.6",alpha=0.55,label="BEFORE"); a_.hist(x1,xb,histtype="step",color="#1f3a93",lw=2.0,label="AFTER")
    if lg: a_.set_xscale("log")
    a_.set_title(f"{t}: median {np.median(x0):.2f} → {np.median(x1):.2f}",fontsize=EB.FT,loc="left"); a_.set_xlabel(xl,fontsize=EB.FS); a_.set_ylabel("count",fontsize=EB.FS); a_.grid(alpha=0.3); a_.tick_params(labelsize=EB.FK)
fig.suptitle(TITLE+"\nfull Delaunay triangulation of the slab |z| < 1.5 Gpc",fontsize=EB.FT+1,fontweight="bold")
fig.legend(handles=[Rectangle((0,0),1,1,fc="0.6",alpha=0.55,label="BEFORE (random start)"),Line2D([],[],color="#1f3a93",lw=2,label="AFTER edits")],loc="lower center",ncol=2,frameon=False,fontsize=EB.FL,bbox_to_anchor=(0.5,0.02))
fig.savefig(f"{OD}/comp_web_hist.png",dpi=130,bbox_inches="tight"); plt.close(fig)
json.dump(info,open(f"{OD}/comp_web_hist.json","w"),indent=1); print(json.dumps(info,indent=1))

"""Pair energies and net forces, foam s2 web (full 3D, inside 30 Gpc), BEFORE vs AFTER the entropy-rule edits.
m = the model's excess mass (final_*.npz 'Mp'): clusters +M; voids m = DELTA_V*rho*V < 0 (DELTA_V = -0.70). U_ij = -G m_i m_j / r_eff,
force on i: F_i = sum_j G m_i m_j (x_j - x_i)/(|x_j-x_i| r_eff^2)  (m_i m_j > 0 attract: cc, vv; < 0 repel: cv).
Contact softening: r_eff = max(r_ij, R_i + R_j) (body radii: voids 1.5 Gpc, clusters 0.158 Gpc).
Pair sets: 'del' = 3D Delaunay neighbours (primary); 'r6' = all pairs within 6 Gpc. Force stats use 'interior' bodies (<= 26 Gpc from us) to avoid the one-sided edge."""
import json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import Delaunay,cKDTree
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle,Rectangle
from matplotlib.lines import Line2D
from matplotlib.colors import SymLogNorm
import edit_before_after as EB, seq_core
from common import *
G=6.674e-11; GP=seq_core.GPC
d=json.load(open("minedit_ent/b_foam_rim_s2_tidalcore_entropy.json")); base=d["tag"].replace("_entropy","")
LAY={"before":EB.layout(base),"after":EB.layout(base,d)}; n=len(LAY["before"][0])
idx=np.array(d["idx"],int); MV=np.zeros(n,bool); MV[idx[np.array(d["moved"],bool)]]=True
def ringdir(X):
    rho=np.hypot(X[:,0],X[:,1]); p=np.c_[X[:,0]/rho*RMIN,X[:,1]/rho*RMIN,np.zeros(len(X))]; u=p-X; return u/np.maximum(np.linalg.norm(u,axis=1),1e-9)[:,None]
R={}
for vw,(X,r,v,M,_) in LAY.items():
    al=np.where((r>0)&(np.linalg.norm(X,axis=1)<=30))[0]; Xa=X[al]; ra=r[al]; va=v[al]; ma=M[al]
    S=Delaunay(Xa).simplices; Ed=np.unique(np.sort(np.vstack([S[:,[a,b]] for a in range(4) for b in range(a+1,4)]),axis=1),axis=0)
    Er=np.array(sorted(cKDTree(Xa).query_pairs(6.0)))
    out={}
    for nm,E in (("del",Ed),("r6",Er)):
        dv=Xa[E[:,1]]-Xa[E[:,0]]; rr=np.linalg.norm(dv,axis=1); re=np.maximum(rr,ra[E[:,0]]+ra[E[:,1]])
        mm=ma[E[:,0]]*ma[E[:,1]]; U=-G*mm/(re*GP); f=(G*mm/(re*GP)**2)[:,None]*dv/rr[:,None]      # force on E0 from E1
        F=np.zeros((len(al),3)); np.add.at(F,E[:,0],f); np.add.at(F,E[:,1],-f)
        ty=np.where(va[E[:,0]]&va[E[:,1]],"vv",np.where(~va[E[:,0]]&~va[E[:,1]],"cc","cv"))
        out[nm]=dict(E=E,U=U,ty=ty,F=F)
    R[vw]=dict(al=al,X=Xa,v=va,r=ra,out=out)
info={}
for vw in R:
    for nm in ("del","r6"):
        o=R[vw]["out"][nm]; info.setdefault(nm,{})[vw]={t:dict(n=int((o["ty"]==t).sum()),total_J=float(o["U"][o["ty"]==t].sum()),median_J=float(np.median(o["U"][o["ty"]==t]))) for t in ("cc","cv","vv")}
        info[nm][vw]["total_J"]=float(o["U"].sum())
        X=R[vw]["X"]; F=o["F"]; inn=np.linalg.norm(X,axis=1)<=26; mv=MV[R[vw]["al"]]; Fr=np.sum(F*X/np.linalg.norm(X,axis=1)[:,None],1); Fg=np.sum(F*ringdir(X),1); Fm=np.linalg.norm(F,axis=1)
        o.update(Fm=Fm,Fr=Fr,Fg=Fg,inn=inn,mv=mv)
        st={}
        for lab,m in (("clusters",~R[vw]["v"]),("voids",R[vw]["v"]),("moved",mv),("not moved",~mv)):
            m=m&inn; st[lab]=dict(n=int(m.sum()),median_absF_N=float(np.median(Fm[m])),mean_F_radial_N=float(np.mean(Fr[m])),mean_F_toward_ring_N=float(np.mean(Fg[m])),
                                  frac_toward_ring=float(np.mean(Fg[m]>0)))
        near=inn&(np.abs(np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2<9)
        st["within 3 Gpc of ring"]=dict(n=int(near.sum()),median_absF_N=float(np.median(Fm[near])),mean_F_toward_ring_N=float(np.mean(Fg[near])),frac_toward_ring=float(np.mean(Fg[near]>0)))
        info[nm][vw]["force"]=st
# moved bodies: same bodies before vs after (Delaunay)
def per_body(vw,nm):
    o=R[vw]["out"][nm]; full=np.full(n,np.nan); a=R[vw]["al"]; return {k:(lambda arr: (full.copy().__setitem__(a,arr) or None))(o[k]) for k in ()}
json.dump(info,open("pair_energy.json","w"),indent=1)
# ---------- (1) energy histograms
cols={"before":"0.55","after":"#1f3a93"}; TY=[("cc","cluster–cluster (attract)"),("cv","cluster–void (repel)"),("vv","void–void (attract)")]
fig,ax=plt.subplots(3,1,figsize=(8.5,13)); fig.subplots_adjust(hspace=0.75)
for a_,(t,lab) in zip(ax,TY):
    allU=np.concatenate([R[vw]["out"]["del"]["U"][R[vw]["out"]["del"]["ty"]==t] for vw in R]); lo,hi=np.percentile(np.abs(allU),[0.2,99.9])
    bins=np.logspace(np.log10(lo),np.log10(hi),45)
    for vw in R:
        o=R[vw]["out"]["del"]; U=o["U"][o["ty"]==t]; s=info["del"][vw][t]
        a_.hist(np.abs(U),bins,histtype="step",color=cols[vw],lw=2.2,label=f"{vw.upper()}: {s['n']:,} pairs, total {s['total_J']:.2e} J, median {s['median_J']:.2e} J")
    a_.set_xscale("log"); a_.set_xlabel(f"|U_ij| (J); sign: U {'< 0' if t!='cv' else '> 0'} for every {t} pair",fontsize=12); a_.set_ylabel("pairs",fontsize=12)
    a_.set_title(f"{lab}, Delaunay neighbours",fontsize=13); a_.tick_params(labelsize=11); a_.grid(alpha=0.3)
    a_.legend(loc="upper center",bbox_to_anchor=(0.5,-0.22),frameon=False,fontsize=10.5)
fig.suptitle("Pair energies U = −G m_i m_j / r (contact-softened), foam s2",fontsize=13.5,fontweight="bold",y=0.93)
fig.savefig("energy_hist_by_type.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# ---------- (2) force histograms
fig,ax=plt.subplots(3,1,figsize=(8.5,13)); fig.subplots_adjust(hspace=0.75)
for a_,(lab,key) in zip(ax[:2],(("clusters",False),("voids",True))):
    allF=np.concatenate([R[vw]["out"]["del"]["Fm"][(R[vw]["v"]==key)&R[vw]["out"]["del"]["inn"]] for vw in R]); bins=np.logspace(np.log10(np.percentile(allF,0.5)),np.log10(allF.max()),45)
    for vw in R:
        o=R[vw]["out"]["del"]; m=(R[vw]["v"]==key)&o["inn"]; a_.hist(o["Fm"][m],bins,histtype="step",color=cols[vw],lw=2.2,label=f"{vw.upper()}: median |F| {np.median(o['Fm'][m]):.2e} N ({m.sum():,} bodies)")
    a_.set_xscale("log"); a_.set_xlabel("|net force| (N)",fontsize=12); a_.set_ylabel("bodies",fontsize=12); a_.set_title(f"Net force on {lab} (Delaunay neighbours, ≤ 26 Gpc from us)",fontsize=13)
    a_.tick_params(labelsize=11); a_.grid(alpha=0.3); a_.legend(loc="upper center",bbox_to_anchor=(0.5,-0.22),frameon=False,fontsize=11)
a_=ax[2]
allg=np.concatenate([R[vw]["out"]["del"]["Fg"][R[vw]["out"]["del"]["inn"]&R[vw]["out"]["del"]["mv"]] for vw in R]); L_=np.percentile(np.abs(allg),99); bins=np.linspace(-L_,L_,41)
for vw in R:
    o=R[vw]["out"]["del"]; m=o["inn"]&o["mv"]; g=o["Fg"][m]
    a_.hist(np.clip(g,-L_,L_),bins,histtype="step",color=cols[vw],lw=2.2,label=f"{vw.upper()}: {100*np.mean(g>0):.0f}% pushed toward the ring, mean {g.mean():+.2e} N")
a_.axvline(0,color="0.3",lw=1); a_.set_xlabel("force component toward the nearest point of the CMB ring (N); > 0 = toward the ring",fontsize=11.5); a_.set_ylabel("bodies",fontsize=12)
a_.set_title(f"Moved bodies only ({int((R['after']['out']['del']['inn']&R['after']['out']['del']['mv']).sum())}): pushed toward or away from the ring?",fontsize=13); a_.tick_params(labelsize=11); a_.grid(alpha=0.3)
a_.legend(loc="upper center",bbox_to_anchor=(0.5,-0.22),frameon=False,fontsize=11)
fig.suptitle("Net forces F_i = Σ_j pair forces, foam s2",fontsize=13.5,fontweight="bold",y=0.93)
fig.savefig("net_force_hist.png",dpi=140,bbox_inches="tight"); plt.close(fig)
# ---------- (3) map coloured by force toward ring
g_all=np.concatenate([R[vw]["out"]["del"]["Fg"] for vw in R]); lim=np.percentile(np.abs(g_all),97); lt=np.percentile(np.abs(g_all),20); norm=SymLogNorm(linthresh=lt,vmin=-lim,vmax=lim); cmap=plt.get_cmap("RdBu_r")
fig=plt.figure(figsize=(13,11.4)); GS=fig.add_gridspec(3,2,height_ratios=[1.0,0.035,0.42],hspace=0.42,wspace=0.1,left=0.06,right=0.98,top=0.9,bottom=0.07)
for k,vw in enumerate(("before","after")):
    a_=fig.add_subplot(GS[0,k]); o=R[vw]["out"]["del"]; X=R[vw]["X"]; v=R[vw]["v"]; r=R[vw]["r"]; sl=np.abs(X[:,2])<1.5; c=cmap(norm(o["Fg"]))
    cp=Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="none",zorder=0); a_.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=a_.transData))
    m=sl&v; a_.add_collection(PatchCollection([Circle(x[:2],rr) for x,rr in zip(X[m],r[m])],fc=c[m],ec=np.where(o["mv"][m][:,None],np.array([[0.04,0.33,0.27,1]]),np.array([[0.5,0.5,0.5,0.6]])),lw=np.where(o["mv"][m],1.6,0.3),zorder=2))
    m=sl&~v; a_.scatter(X[m,0],X[m,1],s=np.where(o["mv"][m],26,12),c=c[m],edgecolor=np.where(o["mv"][m],"#0b5345","0.4"),lw=np.where(o["mv"][m],1.0,0.3),zorder=3)
    a_.add_patch(Circle((0,0),RMIN,fill=False,ec="k",lw=2.0,zorder=6)); a_.plot(0,0,"k+",ms=11,mew=2,zorder=9); frame30(a_,lw=1.2,fs=EB.FK)
    a_.set_xticks([-30,0,30]); a_.set_yticks([-30,0,30]); a_.tick_params(labelsize=EB.FK)
    st=info["del"][vw]["force"]["moved"]; a_.set_title(f"{vw.upper()}: moved bodies {100*st['frac_toward_ring']:.0f}% pushed toward the ring",fontsize=EB.FT)
sm=plt.cm.ScalarMappable(norm=norm,cmap=cmap); cb=fig.colorbar(sm,cax=fig.add_subplot(GS[1,:]),orientation="horizontal"); cb.set_label("net force component toward the CMB ring (N): red = pushed toward the ring, blue = pushed away, white ≈ balanced",fontsize=EB.FK); cb.ax.tick_params(labelsize=EB.FK-1)
EB.residual_below(fig,GS[2,:],d,LAY["before"][4])
fig.legend(handles=[Line2D([],[],marker="o",ls="",ms=6,mfc="w",mec="0.4",label="cluster (dot) / void (disc), slab |z| < 1.5 Gpc"),Line2D([],[],marker="o",ls="",ms=7,mfc="w",mec="#0b5345",mew=1.6,label="dark-green outline = moved body"),
    Line2D([],[],color="k",lw=2,label="CMB ring / our sphere"),Rectangle((0,0),1,1,fc="#fadbd8",label="CoE core cap"),Line2D([],[],color="#a93226",lw=2,marker=">",label="vector from us toward the CoE"),
    Rectangle((0,0),1,1,fc=EB.BAND,ec="#58d68d",ls="--",lw=1,label="±1× CMB band (pale mint)")],loc="upper center",ncol=3,frameon=False,fontsize=EB.FL,bbox_to_anchor=(0.5,0.03))
fig.suptitle("foam s2 — net pair force toward the CMB ring (Delaunay neighbours), BEFORE | AFTER",fontsize=EB.FT+1,fontweight="bold",color="#8e44ad")
fig.savefig("net_force_map.png",dpi=140,bbox_inches="tight")
for nm in ("del","r6"):
    for vw in ("before","after"):
        I=info[nm][vw]; print(nm,vw,"total %.3e"%I["total_J"],{t:(I[t]["n"],"%.2e"%I[t]["total_J"],"%.2e"%I[t]["median_J"]) for t in ("cc","cv","vv")})
        for lab,s in I["force"].items(): print("    ",lab,{a:(round(b,3) if abs(b)<1e3 and isinstance(b,float) else ("%.2e"%b if isinstance(b,float) else b)) for a,b in s.items()})

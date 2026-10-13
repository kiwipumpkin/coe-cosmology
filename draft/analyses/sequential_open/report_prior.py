"""Report for prior_add.py fits (uniform-in-volume prior + closeness penalty): table, distance-from-ring CDFs, leftover vs log-prob, side views."""
import os,json,glob,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from scipy.stats import gamma
from common import *
plt.rcParams.update({"font.size":7.5})
BG=[("none_s0","no background (SR + central mass)","black"),("shallow_cmb_aware_s2_tidalcore","CMB-aware shallow web","#b48fc9"),("b_random_s2_tidalcore","random web","#7f8c8d"),
    ("b_foam_rim_s3_tidalcore","foam web","#5dade2"),("wave_s1","radial wave 150 Mpc","#d35400"),("lattice_s1","quasi-lattice 120 Mpc","#16a085")]
LAMS=[(0,":","unpenalized (λ=0), uniform bank"),(1,"-","prior + penalty, λ=1 (Bayesian)"),(20,"--","λ=20"),(300,(0,(5,1,1,1)),"λ=300"),(3000,(0,(1,1)),"λ=3000")]
def dring(P): P=np.atleast_2d(P); return np.sqrt((np.hypot(P[:,0],P[:,1])-RMIN)**2+P[:,2]**2)
rng=np.random.default_rng(1); U=rng.uniform(-30,30,(4_000_000,3)); s=np.linalg.norm(U,axis=1); U=U[(s<=30)&(s>=0.5)]; Uq=np.sort(dring(U))
Fd=lambda d: np.clip(np.searchsorted(Uq,d)/len(Uq),1/len(Uq),1)
D={}
for t,_,_ in BG:
    for lam,_,_ in LAMS:
        f=f"{OUT}/web_add/{t}_prior_lam{lam}.json"
        if os.path.exists(f): D[(t,lam)]=json.load(open(f))
    f=f"{OUT}/web_add/{t}_omp0_ex1.5.json"     # earlier unpenalized fit with the ring-weighted bank (bodies >= 1.5 Gpc off the ring)
    if os.path.exists(f):
        e=json.load(open(f)); X=np.array(e["X"]); dX=dring(X); lp=float(np.log(Fd(dX)).sum()); N=len(X)
        D[(t,"ex")]=dict(left=e["left"],N=N,n_mass=e["n_mass"],n_void=e["n_void"],within_3gpc=int((dX<3).sum()),within_3gpc_expected=float(N*Fd(np.array([3.0]))[0]),
                         logp_closeness=lp,logp_expected=-N,p_value=float(gamma.sf(-lp,N)),dist_ring=dX.tolist(),final=e["final"],X=e["X"],M=e["M"],resid=e["resid"],stop="earlier fit")
rows=[]
for (t,lam),d in D.items():
    rows.append(dict(background=t,lam=lam,left=round(d["left"],2),N=d["N"],masses=d["n_mass"],voids=d["n_void"],within_3gpc=d["within_3gpc"],within_3gpc_expected=round(d["within_3gpc_expected"],1),
                     logp=round(d["logp_closeness"],1),logp_expected=d["logp_expected"],log10_p=round(np.log10(max(d["p_value"],1e-300)),1),m13=round(d["final"]["rms_13up"],2),m1_12=round(np.sqrt(d["final"]["rms_m1_5"]**2+d["final"]["rms_m6_12"]**2),2),sym=round(d["final"]["sym"],2),stop=d["stop"]))
json.dump(rows,open(f"{OUT}/web_add/prior_table.json","w"),indent=1)
for r in sorted(rows,key=lambda r:(r["background"],str(r["lam"]))): print(r)
# CDF of distance from the ring vs uniform-volume expectation
fig,axs=plt.subplots(3,2,figsize=(7.2,8.6),sharex=True,sharey=True); dd=np.logspace(-1,np.log10(45),300)
for ax,(t,nm,col) in zip(axs.ravel(),BG):
    ax.plot(dd,Fd(dd),color="green",lw=2,label="uniform in volume (prior)")
    for lam,ls,lab in LAMS:
        if (t,lam) in D: x=np.sort(D[(t,lam)]["dist_ring"]); ax.step(x,np.arange(1,len(x)+1)/len(x),where="post",color=col,ls=ls,lw=1.3)
    if (t,"ex") in D: x=np.sort(D[(t,"ex")]["dist_ring"]); ax.step(x,np.arange(1,len(x)+1)/len(x),where="post",color="0.5",lw=0.9,ls="-.")
    ax.axvline(3,color="0.6",lw=0.6,ls=":"); ax.set_xscale("log"); ax.set_title(nm,color=col,fontsize=8,fontweight="bold"); ax.set_ylim(0,1.02)
for ax in axs[-1]: ax.set_xlabel("distance of body centre from the CMB ring (Gpc)")
for ax in axs[:,0]: ax.set_ylabel("fraction of added bodies\ncloser than this")
hd=[Line2D([],[],color="green",lw=2,label="expected if placed uniformly in volume")]+[Line2D([],[],color="0.2",ls=ls,lw=1.3,label=lab) for _,ls,lab in LAMS]+[Line2D([],[],color="0.5",ls="-.",label="earlier fit (ring-weighted bank, ≥1.5 Gpc)"),Line2D([],[],color="0.6",ls=":",lw=0.6,label="3 Gpc")]
fig.legend(handles=hd,loc="lower center",ncol=2,frameon=False,fontsize=7); fig.tight_layout(rect=(0,0.08,1,1)); fig.savefig(f"{OUT}/prior_distance_cdf.png",dpi=200,bbox_inches="tight"); plt.close(fig)
# leftover vs closeness log-probability
fig,ax=plt.subplots(figsize=(7,4.6)); MK={0:"o",1:"s",20:"^",300:"v",3000:"P","ex":"D"}
for (t,lam),d in D.items():
    col=[c for tt,_,c in BG if tt==t][0]; ax.scatter(-d["logp_closeness"]/max(d["N"],1),d["left"],c=col,marker=MK[lam],s=36,edgecolor="k",lw=0.4)
ax.axhline(1,color="green",lw=1.5); ax.axvline(1,color="0.5",ls=":",lw=0.8); ax.text(1.02,0.6,"uniform-volume\nexpectation (1 per body)",fontsize=6.5,color="0.4")
ax.set_xlabel("closeness penalty per body, −log F(d) averaged (prior expects 1.0)"); ax.set_ylabel("final leftover (× 1×)"); ax.set_yscale("log"); ax.set_xscale("log")
hd=[Line2D([],[],marker="o",ls="",color=c,mec="k",label=n) for _,n,c in BG]+[Line2D([],[],marker=MK[k],ls="",mfc="white",mec="k",label=l) for k,l in ((0,"λ=0"),(1,"λ=1"),(20,"λ=20"),(300,"λ=300"),(3000,"λ=3000"),("ex","earlier ring-weighted fit"))]+[Line2D([],[],color="green",lw=1.5,label="1× CMB level")]
fig.legend(handles=hd,loc="lower center",ncol=4,frameon=False,fontsize=6.8); fig.tight_layout(rect=(0,0.16,1,1)); fig.savefig(f"{OUT}/prior_left_vs_logp.png",dpi=200,bbox_inches="tight"); plt.close(fig)
# side views (lambda=1), modeled region only
fig=plt.figure(figsize=(10.5,15.5)); G=fig.add_gridspec(3,2,hspace=0.42,wspace=0.25)
for k,(t,nm,col) in enumerate(BG):
    d=D.get((t,1)); sub=G[k//2,k%2].subgridspec(2,1,height_ratios=[1.7,0.8],hspace=0.45); ax=fig.add_subplot(sub[0])
    if d is None: ax.set_title(nm+": running"); ax.axis("off"); continue
    f=np.load(f"{OUT}/final_{t}.npz")
    if "Xp" in f: v=f["kp"]=="void"; Xv,rv,Xc=f["Xp"][v],f["rp"][v],f["Xp"][~v]
    elif "Rv" in f: Xv=f["X"]; rv=np.full(len(Xv),float(f["Rv"])); Xc=np.zeros((0,3))
    else: Xv=np.zeros((0,3)); rv=np.zeros(0); Xc=np.zeros((0,3))
    sv=(np.abs(Xv[:,2])<1.5)&(np.linalg.norm(Xv,axis=1)<=30); sc=(np.abs(Xc[:,2])<1.5)&(np.linalg.norm(Xc,axis=1)<=30)
    cp=Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="none",zorder=0); ax.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=ax.transData))
    ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(Xv[sv],rv[sv])],fc="#d6eaf8",ec="#85c1e9",lw=0.15,zorder=2))
    if len(Xc): ax.scatter(Xc[sc,0],Xc[sc,1],s=0.5,c="#e6b0aa",lw=0,zorder=3)
    if "zoom" in f: ax.text(0,-29,"background field fills the region (see fig_wave_lattice_alone.png)",fontsize=5.5,ha="center",color=col)
    X=np.array(d["X"]); M=np.array(d["M"]); vm=M<0; rho=seq_core.rho_gpc(X)*seq_core.GPC**3; Rvd=(np.abs(M)/(0.7*rho*4/3*np.pi))**(1/3)
    ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(X[vm],Rvd[vm])],fc="none",ec="#1a5276",lw=1.0,zorder=6))
    ax.scatter(X[~vm,0],X[~vm,1],s=12+30*np.abs(M[~vm])/2e52,c="#c0392b",edgecolor="k",lw=0.3,zorder=7)
    ax.add_patch(Circle((0,0),RMIN,fill=False,ec="#1f4e79",lw=1.3,zorder=5)); ax.plot(0,0,"k+"); frame30(ax); ax.tick_params(labelsize=5.5)
    for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(2)
    ax.set_title(f"{nm} — prior + penalty (λ=1)\nadded {d['n_mass']} masses + {d['n_void']} voids; within 3 Gpc of ring: {d['within_3gpc']} (uniform expects {d['within_3gpc_expected']:.1f})\nlog-prob of closeness {d['logp_closeness']:.0f} (prior expects {d['logp_expected']}); p = {d['p_value']:.1g}",fontsize=6,color=col,fontweight="bold")
    r_=fig.add_subplot(sub[1]); res=np.array(d["resid"])/FLOOR
    r_.axhspan(-1,1,color="green",alpha=0.3,lw=0); r_.axhline(0,color="0.6",lw=0.5); r_.plot(np.arange(360),res,color=col,lw=0.9)
    lim=1.3*max(np.abs(res).max(),2); r_.set_yscale("symlog",linthresh=1,linscale=0.6); r_.set_ylim(-lim,lim); r_.set_xlim(0,360); r_.set_xticks([0,90,180,270,360]); r_.tick_params(labelsize=5.5)
    fs=d["final"]; r_.set_title(f"left {d['left']:.2f}×  (m1–5 {fs['rms_m1_5']:.2f}, m6–12 {fs['rms_m6_12']:.2f}, 13+ {fs['rms_13up']:.2f}; symmetric {fs['sym']*100:.0f}%)",fontsize=6.3,color=col)
    r_.set_xlabel("ℓ (°) around the ring",fontsize=5.8); r_.set_ylabel("× $10^{-5}$ (symlog)",fontsize=5.5)
hd=[Line2D([],[],marker="o",ls="",mfc="#c0392b",mec="k",ms=5,label="added mass (size ∝ mass; all z, projected)"),Line2D([],[],marker="o",ls="",mfc="none",mec="#1a5276",ms=7,label="added void (to scale)"),
    Line2D([],[],marker="o",ls="",mfc="#d6eaf8",mec="#85c1e9",ms=6,label="background void (slab |z|<1.5)"),Line2D([],[],marker="o",ls="",color="#e6b0aa",ms=3,label="background cluster"),
    Line2D([],[],color="#1f4e79",lw=1.3,label="CMB ring / our sphere"),plt.Rectangle((0,0),1,1,fc="#fadbd8",label="core cap inside 30 Gpc (filled)"),Line2D([],[],color="0.4",ls=":",label="30 Gpc"),
    Line2D([],[],color="#a93226",marker=">",label="direction to the CoE"),plt.Rectangle((0,0),1,1,fc="green",alpha=0.3,label="±1× (CMB level)")]
fig.legend(handles=hd,loc="lower center",ncol=3,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.05)); fig.savefig(f"{OUT}/prior_side_views.png",dpi=170,bbox_inches="tight"); plt.close(fig)

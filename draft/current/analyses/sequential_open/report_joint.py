"""Report: joint multi-start bodies-on-top (joint_add.py) -> joint_*.png, joint_add/summary.json"""
import os,json,glob,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
from common import *
plt.rcParams.update({"font.size":7.5})
BG=[("none_s0","no background (SR + central mass)","black"),("shallow_cmb_aware_s2_tidalcore","CMB-aware shallow web","#b48fc9"),("b_random_s2_tidalcore","random web","#7f8c8d"),
    ("b_foam_rim_s3_tidalcore","foam web","#5dade2"),("wave_s1","radial wave 150 Mpc","#d35400"),("lattice_s1","quasi-lattice 120 Mpc","#16a085")]
def dring(X): X=np.atleast_2d(X); return np.sqrt((np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2)
R_=[json.load(open(f)) for f in glob.glob(f"{OUT}/joint_add/raw/*.json")]
B={}; NST={}
for r in R_:
    k=(r["tag"],r["N"]); NST[k]=NST.get(k,0)+1
    if r["ok"] and (k not in B or r["left"]<B[k]["left"]): B[k]=r
bgs=[b for b in BG if any(k[0]==b[0] for k in B)]
# uniform-in-volume expectation over the allowed region (ball 30, >= 0.5 from us, surface >= 3 Gpc from ring for a point body)
rng=np.random.default_rng(0); U=rng.uniform(-30,30,(2_000_000,3)); s=np.linalg.norm(U,axis=1); U=U[(s<=30)&(s>=0.5)]; U=U[dring(U)>=3]
# earlier fits (solution_space/asym, no-mirror, before any ring rule) for comparison
E=json.load(open("/workspace/clump_scan/solution_space/asym/asym_results.json"))["best"]; EX=np.concatenate([np.array(v["X"]) for v in E.values()])
# 1: best leftover vs N; split
fig,axs=plt.subplots(3,1,figsize=(7,10),sharex=True)
for t,nm,col in bgs:
    Ns=sorted(N for (tt,N) in B if tt==t); b=[B[(t,N)] for N in Ns]
    axs[0].plot(Ns,[x["left"] for x in b],"-o",color=col,ms=4,label=nm)
    axs[1].plot(Ns,[np.hypot(x["spec"]["rms_m1_5"],x["spec"]["rms_m6_12"]) for x in b],"-o",color=col,ms=4); axs[1].plot(Ns,[x["spec"]["rms_13up"] for x in b],"--s",color=col,ms=3.5,mfc="white")
    axs[2].plot(Ns,[x["spec"]["sym"] for x in b],"-o",color=col,ms=4)
for ax in axs[:2]: ax.axhline(1,color="green",lw=1.5); ax.set_yscale("log")
axs[0].set_ylabel("best leftover (× 1×)"); axs[0].set_title("Best fit over random starts (all bodies fitted jointly; surface ≥ 3 Gpc from the ring)",fontsize=8)
axs[1].set_ylabel("part of the leftover (× 1×)"); axs[1].set_title("solid: modes m = 1–12    dashed: m ≥ 13",fontsize=8)
axs[2].set_ylabel("symmetric share of m=1–12"); axs[2].axhline(0.5,color="0.5",ls=":"); axs[2].set_ylim(0,1); axs[2].set_xlabel("number of added bodies N (half masses, half voids)")
hd=[Line2D([],[],color=c,marker="o",label=n) for _,n,c in bgs]+[Line2D([],[],color="green",lw=1.5,label="1× (CMB level)")]
fig.legend(handles=hd,loc="lower center",ncol=3,frameon=False,fontsize=7); fig.tight_layout(rect=(0,0.06,1,1)); fig.savefig(f"{OUT}/joint_best_vs_N.png",dpi=200,bbox_inches="tight"); plt.close(fig)
# 2: distance distributions (best fits, all N pooled), vs uniform-volume and the earlier fits
fig,axs=plt.subplots(2,1,figsize=(7,7.5))
for ax,key,ref,eref,xl in ((axs[0],"dist_us",np.linalg.norm(U,axis=1),np.linalg.norm(EX,axis=1),"distance of body from us (Gpc)"),(axs[1],"dist_ring",dring(U),dring(EX),"distance of body centre from the CMB ring (Gpc)")):
    bins=np.linspace(0,45,46)
    ax.hist(ref,bins=bins,density=True,histtype="stepfilled",color="green",alpha=0.2,label="uniform in volume over the allowed region")
    ax.hist(eref,bins=bins,density=True,histtype="step",color="0.4",lw=1.4,ls="--",label="earlier joint fits (solution_space/asym, best per N)")
    for t,nm,col in bgs:
        v=np.concatenate([B[k][key] for k in B if k[0]==t]); ax.hist(v,bins=bins,density=True,histtype="step",color=col,lw=1.1)
    ax.set_xlabel(xl); ax.set_ylabel("fraction of bodies per Gpc")
axs[1].axvline(3,color="red",lw=0.8,ls=":"); axs[1].text(3.3,axs[1].get_ylim()[1]*0.9,"3 Gpc",color="red",fontsize=6.5)
hd=[plt.Rectangle((0,0),1,1,fc="green",alpha=0.2,label="uniform in volume (allowed region)"),Line2D([],[],color="0.4",ls="--",lw=1.4,label="earlier joint fits (no ring rule, reach 43 Gpc)")]+[Line2D([],[],color=c,label=n) for _,n,c in bgs]
fig.legend(handles=hd,loc="lower center",ncol=2,frameon=False,fontsize=7); fig.tight_layout(rect=(0,0.1,1,1)); fig.savefig(f"{OUT}/joint_distances.png",dpi=200,bbox_inches="tight"); plt.close(fig)
# 3: side views of the best fit at the largest N done
fig=plt.figure(figsize=(10.5,5.3*int(np.ceil(len(bgs)/2)))); G=fig.add_gridspec(int(np.ceil(len(bgs)/2)),2,hspace=0.42,wspace=0.25)
for q,(t,nm,col) in enumerate(bgs):
    Nb=min((N for (tt,N) in B if tt==t),key=lambda N:B[(t,N)]["left"]); d=B[(t,Nb)]
    sub=G[q//2,q%2].subgridspec(2,1,height_ratios=[1.7,0.8],hspace=0.45); ax=fig.add_subplot(sub[0]); f=np.load(f"{OUT}/final_{t}.npz")
    if "Xp" in f: v=f["kp"]=="void"; Xv,rv,Xc=f["Xp"][v],f["rp"][v],f["Xp"][~v]
    elif "Rv" in f: Xv=f["X"]; rv=np.full(len(Xv),float(f["Rv"])); Xc=np.zeros((0,3))
    else: Xv=np.zeros((0,3)); rv=np.zeros(0); Xc=np.zeros((0,3))
    sv=(np.abs(Xv[:,2])<1.5)&(np.linalg.norm(Xv,axis=1)<=30); sc=(np.abs(Xc[:,2])<1.5)&(np.linalg.norm(Xc,axis=1)<=30)
    cp=Circle(C_G[:2],A_IN_G,fc="#fadbd8",ec="none",zorder=0); ax.add_patch(cp); cp.set_clip_path(Circle((0,0),30.0,transform=ax.transData))
    th=np.linspace(0,2*np.pi,300)
    for rr in (RMIN-3,RMIN+3): ax.plot(rr*np.cos(th),rr*np.sin(th),color="red",lw=0.5,ls=":",zorder=4)
    ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(Xv[sv],rv[sv])],fc="#d6eaf8",ec="#85c1e9",lw=0.15,zorder=2))
    if len(Xc): ax.scatter(Xc[sc,0],Xc[sc,1],s=0.5,c="#e6b0aa",lw=0,zorder=3)
    if "zoom" in f: ax.text(0,-29,"background field fills the region (see fig_wave_lattice_alone.png)",fontsize=5.5,ha="center",color=col)
    X=np.array(d["X"]); M=np.array(d["M"]); Rv=np.array(d["Rv"]); vm=M<0
    ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(X[vm],Rv[vm])],fc="#1a5276",alpha=0.25,ec="#1a5276",lw=1.0,zorder=6))
    ax.scatter(X[~vm,0],X[~vm,1],s=14+40*np.abs(M[~vm])/2e52,c="#c0392b",edgecolor="k",lw=0.4,zorder=7)
    ax.add_patch(Circle((0,0),RMIN,fill=False,ec="#1f4e79",lw=1.3,zorder=5)); ax.plot(0,0,"k+"); frame30(ax); ax.tick_params(labelsize=5.5)
    for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(2)
    ax.set_title(f"{nm}: best fit, N = {Nb} ({(~vm).sum()} masses + {vm.sum()} voids, all z projected)\nbackground alone + SR: {d['background']['total']:.1f}×; body surfaces ≥ {min(d['surface_ring']):.1f} Gpc from the ring",fontsize=6.3,color=col,fontweight="bold")
    r_=fig.add_subplot(sub[1]); res=np.array(d["resid"])/FLOOR
    r_.axhspan(-1,1,color="green",alpha=0.3,lw=0); r_.axhline(0,color="0.6",lw=0.5); r_.plot(np.arange(360),res,color=col,lw=0.9)
    lim=1.3*max(np.abs(res).max(),2); r_.set_yscale("symlog",linthresh=1,linscale=0.6); r_.set_ylim(-lim,lim); r_.set_xlim(0,360); r_.set_xticks([0,90,180,270,360]); r_.tick_params(labelsize=5.5)
    fs=d["spec"]; r_.set_title(f"left {d['left']:.2f}×  (m1–5 {fs['rms_m1_5']:.2f}, m6–12 {fs['rms_m6_12']:.2f}, 13+ {fs['rms_13up']:.2f}; symmetric {fs['sym']*100:.0f}%)",fontsize=6.3,color=col)
    r_.set_xlabel("ℓ (°) around the ring",fontsize=5.8); r_.set_ylabel("× $10^{-5}$ (symlog)",fontsize=5.5)
hd=[Line2D([],[],marker="o",ls="",mfc="#c0392b",mec="k",ms=5,label="added mass (size ∝ mass)"),plt.Rectangle((0,0),1,1,fc="#1a5276",alpha=0.25,ec="#1a5276",label="added void (to scale)"),
    Line2D([],[],marker="o",ls="",mfc="#d6eaf8",mec="#85c1e9",ms=6,label="background void (slab |z|<1.5)"),Line2D([],[],marker="o",ls="",color="#e6b0aa",ms=3,label="background cluster"),
    Line2D([],[],color="#1f4e79",lw=1.3,label="CMB ring / our sphere"),Line2D([],[],color="red",ls=":",lw=0.6,label="±3 Gpc from the ring (in-plane)"),plt.Rectangle((0,0),1,1,fc="#fadbd8",label="core cap inside 30 Gpc (filled)"),
    Line2D([],[],color="0.4",ls=":",label="30 Gpc"),Line2D([],[],color="#a93226",marker=">",label="direction to the CoE"),plt.Rectangle((0,0),1,1,fc="green",alpha=0.3,label="±1× (CMB level)")]
fig.legend(handles=hd,loc="upper center",ncol=4,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.0)); fig.savefig(f"{OUT}/joint_side_views.png",dpi=170,bbox_inches="tight"); plt.close(fig)
S={f"{t}|N{N}":dict(left=round(b["left"],2),starts=NST[(t,N)],ok_starts=sum(1 for r in R_ if r["tag"]==t and r["N"]==N and r["ok"]),m1_12=round(np.hypot(b["spec"]["rms_m1_5"],b["spec"]["rms_m6_12"]),2),m13=round(b["spec"]["rms_13up"],2),sym=round(b["spec"]["sym"],2),
     inside_us=b["inside_us"],in_core=b["in_core"],min_surface_ring=round(min(b["surface_ring"]),2),median_dist_us=round(float(np.median(b["dist_us"])),1),median_dist_ring=round(float(np.median(b["dist_ring"])),1)) for (t,N),b in sorted(B.items())}
json.dump(S,open(f"{OUT}/joint_add/summary.json","w"),indent=1)
for k,v in S.items(): print(k,v)
print("uniform-volume medians: dist_us %.1f  dist_ring %.1f ; earlier fits medians: dist_us %.1f dist_ring %.1f"%(np.median(np.linalg.norm(U,axis=1)),np.median(dring(U)),np.median(np.linalg.norm(EX,axis=1)),np.median(dring(EX))))

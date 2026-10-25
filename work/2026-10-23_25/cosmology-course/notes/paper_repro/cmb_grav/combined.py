import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, json
exec(open("cmb_grav.py").read().split("Ms=np.r_")[0])          # reuse geometry/events/leg (no run)
N="/workspace/cosmology-course/notes/"
b=np.loadtxt(N+"migkas-A-vs-angle-bins.csv",delimiter=",",skiprows=1)
mc=np.loadtxt(N+"catwise-exp-migkas-fit-d14.3.csv",delimiter=",",skiprows=1)
t3=[r for r in json.load(open("../cmb_iso3/cmb_iso3_table.json")) if r["beta"]==0 and r["law"]=="power_law" and r["M"]==0][0]
E=events(0); base=np.array([Trec*e["dop"]/e["stretch"] for e in E])
me=np.array([Menc1(e["rnow"]) for e in E]); pe=-Gc2*me/Rrec; po=-Gc2*Menc1(D)/D
def boosted(Tv,beta):
    g=1/np.sqrt(1-beta**2); v=-mu; vo=(v-beta)/(1-beta*v)
    return np.interp(np.arccos(np.clip(-vo,-1,1)),th,Tv)/(g*(1-beta*v))
def dTT(Tv): m=np.trapezoid(Tv*W,th)/2; x=Tv/m-1; return x,float(np.sqrt(np.trapezoid(x*x*W,th)/2))
cases=[("cmb_iso3 β=0, M=0",base),(f"β=0, best SR boost β_obs={t3['best_boost']:.4f} away from CoE",boosted(base,t3["best_boost"]))]
for M in (5.62e15,5.62e16):
    cases.append((f"cmb_grav β=0, M={M:.2e} M☉ weak field"+(" (FLAG: |Φ|/c²≥0.1)" if M>1e16 else ""),base*(1+M*(pe-po))))
l0=15.2; ell=np.arange(0,361); thl=np.radians(np.abs(((ell-l0+180)%360)-180))
fig=plt.figure(figsize=(12,13)); gs=fig.add_gridspec(3,1,height_ratios=[1,1,0.6]); a0=fig.add_subplot(gs[0]); ax=[a0,fig.add_subplot(gs[1],sharex=a0),fig.add_subplot(gs[2])]
ax[0].errorbar(b[:,0],b[:,1],yerr=[b[:,2],b[:,3]],fmt="o",ms=4,capsize=2,color="k",label="Migkas A/A_all bins (migkas-A-vs-angle-bins.csv)")
ax[0].plot(mc[:,0],mc[:,1],color="tab:orange",label="existing model curve, D=58.2, R=14.3, l0=15.16 (catwise-exp-migkas-fit-d14.3.csv; not refit)")
ax[0].axhline(1,color="gray",lw=.7); ax[0].set_ylabel("A / A_all"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
out={}
for (lab,Tv),col in zip(cases,("k","tab:green","tab:blue","tab:red")):
    x,r=dTT(Tv); out[lab]=r; ax[1].plot(ell,np.interp(thl,th,x),color=col,label=f"{lab}; RMS={r:.4g}")
for y,lab,ls in ((5.75e-2,"paper level 5.75e-2 (start)","--"),(2.37e-3,"paper level 2.37e-3 (after SR+mass)",":")):
    for s in (1,-1): ax[1].axhline(s*y,color="purple",ls=ls,lw=1,label=lab if s==1 else None)
ax[1].axvline(l0,color="gray",lw=.6); ax[1].axvline(l0+180,color="gray",lw=.6)
ax[1].set_ylabel("CMB dT/T (full, dipole kept)"); ax[1].set_xlabel("galactic longitude ℓ (deg), in-plane cut through CoE axis ℓ=15.2 (θ=0) / 195.2 (θ=180)")
ax[1].legend(fontsize=7.5,loc="lower right"); ax[1].grid(alpha=.3); ax[1].set_xlim(0,360)
labs=list(out); ax[2].bar(range(len(labs)),[out[k] for k in labs],color=["k","tab:green","tab:blue","tab:red"])
for i,k in enumerate(labs): ax[2].text(i,out[k],f"{out[k]:.4g}",ha="center",va="bottom",fontsize=8)
ax[2].axhline(5.75e-2,color="purple",ls="--",label="paper 5.75e-2 (start)"); ax[2].axhline(2.37e-3,color="purple",ls=":",label="paper 2.37e-3 (after SR+mass)")
ax[2].set_xticks(range(len(labs))); ax[2].set_xticklabels(["β=0, M=0","best SR boost","M=5.62e15","M=5.62e16 (flag)"]); ax[2].set_ylabel("RMS dT/T"); ax[2].legend(fontsize=8); ax[2].grid(alpha=.3,axis="y")
plt.tight_layout(); plt.savefig("combined_aall_cmb.png",dpi=130); print(out)

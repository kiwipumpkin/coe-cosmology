exec(open("paper_repro/stretch/chain_stretch.py").read().split("D,Rv,phi0")[0])
import numpy as np
D,Rv,phi0=58.2,14.3,15.16
rhit=np.sqrt(D**2+Rv**2+2*D*Rv*np.cos(np.radians(ELL-phi0)))
Gs=G_FULL*(rhit/rhit.mean())**(-2/3); Gs/=Gs.mean()
def fit(G,mu0):
    g=lambda p: f([p[0],p[1],0.0,p[2],p[3]],G) if mu0 else f(p,G)
    bd=[(0.01,0.99),(0.01,3),(1,3),(0,360)] if mu0 else [(0.01,0.99),(0.01,3),(0,0.5),(1,3),(0,360)]
    de=differential_evolution(g,bd,seed=1,tol=1e-10,maxiter=800,polish=False)
    r=minimize(g,de.x,method="Nelder-Mead",options=dict(xatol=1e-10,fatol=1e-13,maxiter=40000))
    p=list(r.x); p=[p[0],p[1],0.0,p[2],p[3]] if mu0 else p
    return np.array(p),r.fun
C={}
for nm,G,m0 in (("paper chain (with mass)",G_FULL,0),("paper chain + stretch (with mass)",Gs,0),("paper chain, zero mass",G_FULL,1),("paper chain + stretch, zero mass",Gs,1)):
    p,r=fit(G,m0)
    if nm=="paper chain (with mass)" and f(P2full,G)<r: p,r=np.array(P2full),f(P2full,G)
    b,q,mu,s,ph=p; th=np.radians(ELL-ph); re0=np.sqrt(1+q*q+2*q*np.cos(th))
    vo=np.sqrt(b*b+2*mu*(1/s-1)); ve=np.sqrt(b*b*re0**2+2*mu*(1/(s*re0)-1/re0))
    R=R_of(G,p); C[nm]=R
    print(f"{nm}: RMS={r:.4e} max|R|={np.abs(R).max():.3e} beta1={b:.5f} q={q:.5f} mu={mu:.5f} s={s:.4f} phi={ph%360:.3f} vo={vo:.3f} ve={ve.min():.3f}-{ve.max():.3f} max2mu/r={max(2*mu/s,(2*mu/(s*re0)).max()):.3f}")
b_=np.loadtxt("migkas-A-vs-angle-bins.csv",delimiter=",",skiprows=1)
plt.rcParams.update({"font.size":17}); fig,ax=plt.subplots(2,1,figsize=(13,12),sharex=True)
ax[0].errorbar(b_[:,0],b_[:,1],yerr=[b_[:,2],b_[:,3]],fmt="o",color="k",ms=5,capsize=2,label="Migkas A/A_all bins")
ax[0].plot(ELL,full,lw=3,color="tab:orange",label="paper WLS geometry D=58.2, R=14.3, ℓ0=15.2")
ax[0].set_ylabel("A / A_all"); ax[0].legend(fontsize=13); ax[0].grid(alpha=.3)
for (nm,R),st in zip(C.items(),(dict(color="tab:red",lw=3),dict(color="tab:blue",lw=3,ls="--"),dict(color="tab:green",lw=3),dict(color="tab:purple",lw=3,ls="--"))):
    ax[1].plot(ELL,R*1e3,label=f"{nm}: RMS = {rms(R):.3e}",**st)
ax[1].axhline(0,color="gray"); ax[1].set_ylabel("CMB residual R × 10³"); ax[1].set_xlabel("galactic longitude ℓ (deg)")
lo=min(R.min() for R in C.values())*1e3; hi=max(R.max() for R in C.values())*1e3
ax[1].set_ylim(lo-(hi-lo)*0.45,hi*1.05); ax[1].legend(fontsize=13,loc="lower center",ncol=2); ax[1].grid(alpha=.3); ax[1].set_xlim(0,359)
plt.tight_layout(); plt.savefig("paper_repro/stretch/paper_chain_stretch_mass0.png",dpi=120)

exec(open("paper_repro/stretch/chain_stretch.py").read().split("D,Rv,phi0")[0])
import numpy as np
def S_gr(beta1,q,mu,s,phi,ell=ELL):
    th=np.radians(ell-phi); xe0=1+q*np.cos(th); ye0=q*np.sin(th); re0=np.hypot(xe0,ye0)
    ro,re=s,s*re0
    # exact radial Schwarzschild geodesics: conserved E, local speed vs static observers
    Eo2=(1-2*mu)/(1-beta1**2); Ee2=(1-2*mu/re0)/(1-(beta1*re0)**2)
    vo2=1-(1-2*mu/ro)/Eo2; ve2=1-(1-2*mu/re)/Ee2
    if np.any(np.asarray(1-(beta1*re0)**2)<=0) or vo2<=0 or vo2>=1 or np.any(ve2<=0) or np.any(ve2>=1): return None
    vo,ve=np.sqrt(vo2),np.sqrt(ve2); xe,ye=s*xe0,s*ye0; dx,dy=ro-xe,-ye
    a=vo2-1; b=2*dx*vo; c=dx**2+dy**2; tf=(-b-np.sqrt(b*b-4*a*c))/(2*a)
    nx,ny=(dx+vo*tf)/tf,dy/tf; bex,bey=ve*xe/re,ve*ye/re
    go,ge=1/np.sqrt(1-vo2),1/np.sqrt(1-ve2)
    Dm=go*(1-vo*nx)/(ge*(1-(bex*nx+bey*ny))); Dg=np.sqrt((1-2*mu/re)/(1-2*mu/ro))
    return Dm*Dg,Dm,Dg
p1=np.array(P2full,float); p0=p1.copy(); p0[2]=0
out={}
for nm,p,fn in (("weak-field, mass",p1,S_model),("exact GR, mass",p1,S_gr),("exact GR, zero mass",p0,S_gr)):
    o=fn(*p); out[nm]=resid(G_FULL,o[0]); print(nm,rms(out[nm]))
b_=np.loadtxt("migkas-A-vs-angle-bins.csv",delimiter=",",skiprows=1)
plt.rcParams.update({"font.size":18}); fig,ax=plt.subplots(2,1,figsize=(13,12),sharex=True)
ax[0].errorbar(b_[:,0],b_[:,1],yerr=[b_[:,2],b_[:,3]],fmt="o",color="k",ms=5,capsize=2,label="Migkas bins")
ax[0].plot(ELL,full,lw=3,color="tab:orange",label="paper geometry"); ax[0].set_ylabel("A / A_all"); ax[0].legend(fontsize=14); ax[0].grid(alpha=.3)
for (nm,R),c,ls in zip(out.items(),("tab:red","tab:blue","tab:green"),("-","--","-")):
    ax[1].plot(ELL,R*1e3,color=c,ls=ls,lw=3,label=f"{nm}: RMS {rms(R):.2e}")
ax[1].axhline(0,color="gray"); ax[1].set_ylabel("CMB ΔT/T × 10³"); ax[1].set_xlabel("galactic longitude ℓ (deg)")
ax[1].legend(fontsize=14); ax[1].grid(alpha=.3); ax[1].set_xlim(0,359)
plt.tight_layout(); plt.savefig("paper_repro/stretch/gr_exact.png",dpi=120)

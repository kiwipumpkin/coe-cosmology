import numpy as np, json, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
src=open("/workspace/cosmology-course/notes/cmb-cancel.py").read().split("def obj(")[0]; exec(src)
P2=(0.5866,0.2303,0.09332,3.0,15.2); P2full=(0.5865994328877997,0.23029859560942934,0.09332458226714006,3.0,15.162930413775014)
def R_of(G,p):
    o=S_model(*p); return None if o is None else resid(G,o[0])
def f(p,G):
    if not (0<p[0]<1 and 0<p[1] and p[2]>=0 and 1<=p[3]<=3): return 1.0
    r=R_of(G,p); return 1.0 if r is None else rms(r)
D,Rv,phi0=58.2,14.3,15.16
rhit=np.sqrt(D**2+Rv**2+2*D*Rv*np.cos(np.radians(ELL-phi0)))   # hit radius: θ=0 (ℓ=ϕ) points away from CoE
out={"rms_G":rms(G_FULL-1),"paper_rounded":rms(R_of(G_FULL,P2)),"paper_fullprec":rms(R_of(G_FULL,P2full))}
fits={}
for name,expo in (("stretch_2/3 (fixed-width annulus)",2/3),("stretch_1 (homologous)",1.0)):
    St=(rhit/rhit.mean())**(-expo); G2=G_FULL*St; G2/=G2.mean()
    best=None
    for x0 in (P2full,):
        r=minimize(f,x0,args=(G2,),method="Nelder-Mead",options=dict(xatol=1e-10,fatol=1e-12,maxiter=20000))
        best=r if best is None or r.fun<best.fun else best
    de=differential_evolution(f,[(0.01,0.99),(0.01,3),(0,0.5),(1,3),(0,360)],args=(G2,),seed=1,tol=1e-10,maxiter=600,polish=False)
    if de.fun<best.fun: best=minimize(f,de.x,args=(G2,),method="Nelder-Mead",options=dict(xatol=1e-10,fatol=1e-12,maxiter=20000))
    fits[name]=dict(G=G2,p=best.x,rms=best.fun,rms_G=rms(G2-1),rms_at_paper_params=f(P2full,G2),maxR=float(np.abs(R_of(G2,best.x)).max()))
    print(name,"rms(G-1)=%.4e"%fits[name]["rms_G"],"at paper params %.4e"%fits[name]["rms_at_paper_params"],"refit rms=%.4e max|R|=%.3e"%(best.fun,fits[name]["maxR"]),"params",np.round(best.x,5))
print(out)
Rp=R_of(G_FULL,P2full); k="stretch_2/3 (fixed-width annulus)"; Rs=R_of(fits[k]["G"],fits[k]["p"])
b=np.loadtxt("migkas-A-vs-angle-bins.csv",delimiter=",",skiprows=1)
plt.rcParams.update({"font.size":17}); fig,ax=plt.subplots(2,1,figsize=(13,11),sharex=True)
ax[0].errorbar(b[:,0],b[:,1],yerr=[b[:,2],b[:,3]],fmt="o",color="k",ms=5,capsize=2,label="Migkas A/A_all bins")
ax[0].plot(ELL,full,lw=3,color="tab:orange",label="paper WLS geometry D=58.2, R=14.3, ℓ0=15.2")
ax[0].set_ylabel("A / A_all"); ax[0].legend(fontsize=13); ax[0].grid(alpha=.3)
ax[1].plot(ELL,Rp*1e3,lw=3,color="tab:red",label=f"paper chain: RMS = {rms(Rp):.3e}")
ax[1].plot(ELL,Rs*1e3,lw=3,ls="--",color="tab:blue",label=f"paper chain + shell stretch: RMS = {rms(Rs):.3e}")
ax[1].axhline(0,color="gray"); ax[1].set_ylabel("CMB residual R × 10³"); ax[1].set_xlabel("galactic longitude ℓ (deg)")
ax[1].set_ylim(-6,4.2); ax[1].legend(fontsize=14,loc="lower center",ncol=2); ax[1].grid(alpha=.3); ax[1].set_xlim(0,359)
plt.tight_layout(); plt.savefig("paper_repro/stretch/paper_chain_stretch.png",dpi=120)
json.dump(dict(out,**{n:dict(params=list(map(float,v["p"])),rms=float(v["rms"]),rms_G=v["rms_G"],rms_at_paper_params=v["rms_at_paper_params"],maxR=v["maxR"]) for n,v in fits.items()}),open("paper_repro/stretch/result.json","w"),indent=1)

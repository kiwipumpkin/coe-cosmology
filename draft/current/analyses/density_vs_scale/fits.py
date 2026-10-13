import sys, json, numpy as np
sys.path.insert(0,"/workspace/clump_scan/inside_visible/smooth_field")
from lcdm import sigma_TH
V=lambda R: 4/3*np.pi*np.asarray(R,float)**3
def fit(R,d,tf=lambda d:d):
    x=np.log10(V(R)); y=np.log10(tf(np.asarray(d,float)))
    n=len(x); X=np.c_[np.ones(n),x]; b,res,_,_=np.linalg.lstsq(X,y,rcond=None)
    r=y-X@b; dof=n-2
    s2=(r@r)/dof if dof>0 else np.nan; cov=s2*np.linalg.inv(X.T@X)
    return dict(logA=b[0],alpha=-b[1],se_alpha=np.sqrt(cov[1,1]),rms_dex=np.sqrt(r@r/n),n=n,cov=cov.tolist())
over={"Coma":(2.84,645),"Shapley core":(11.4,20),"Saraswati":(36.6,1.62),"SGW":(79,1.0),"Shapley whole":(43.7,0.8)}
under={"Bootes":(44,0.75),"SDSS voids":(40,0.4),"Eridanus":(314,0.14),"Local Hole":(250,0.24)}
def pts(D,keys): return np.array([D[k][0] for k in keys]),np.array([D[k][1] for k in keys])
F={}
F["over4"]=fit(*pts(over,["Coma","Shapley core","Saraswati","SGW"]))
F["over5"]=fit(*pts(over,list(over)))
F["over_noComa"]=fit(*pts(over,["Shapley core","Saraswati","SGW","Shapley whole"]))
F["over_superclusters"]=fit(*pts(over,["Saraswati","SGW","Shapley whole"]))
F["under_noLH"]=fit(*pts(under,["Bootes","SDSS voids","Eridanus"]))
F["under_LH"]=fit(*pts(under,list(under)))
sat=lambda d:d/(1-d)
F["under_sat_noLH"]=fit(*pts(under,["Bootes","SDSS voids","Eridanus"]),sat)
F["under_sat_LH"]=fit(*pts(under,list(under)),sat)
# LCDM slope dln sigma/dln V
for R in [100,300,1000,2000,4000]:
    e=1.05; s=np.log(sigma_TH(R*e/1000)/sigma_TH(R/1000))/np.log(e**3); F[f"lcdm_slope_{R}"]=-s
def pred(f,R,satur=False):
    x=np.log10(V(R)); y=f["logA"]-f["alpha"]*x
    c=np.array(f["cov"]); sy=np.sqrt(c[0,0]+x*x*c[1,1]+2*x*c[0,1]); v=10**y
    lo,hi=10**(y-sy),10**(y+sy)
    if satur: v,lo,hi=[q/(1+q) for q in (v,lo,hi)]
    return v,lo,hi
peak={1000:0.068,2000:0.055,4000:0.059}; rms={1000:0.0105,2000:0.0124,4000:0.0141}
out={}
for k,f in F.items():
    if not isinstance(f,dict): print(k,round(f,3)); continue
    print(f"{k}: |d| = 10^{f['logA']:.2f} V^-{f['alpha']:.3f} ±{f['se_alpha']:.3f}  (n={f['n']}, scatter {f['rms_dex']:.2f} dex)")
    for R in [1000,2000,4000]:
        v,lo,hi=pred(f,R,k.startswith("under_sat"))
        print(f"   R={R}: {v:.2e} [{lo:.1e},{hi:.1e}]  peak/fit={peak[R]/v:.3g}  rms/fit={rms[R]/v:.3g}  3sig={3*sigma_TH(R/1000):.1e}")
    for R in [830,3750]:
        v,_,_=pred(f,R,k.startswith("under_sat")); print(f"   R={R}: {v:.2e}")
json.dump(F,open("fits.json","w"),indent=1)

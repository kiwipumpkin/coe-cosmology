"""Pure point mass at CoE (no shell, Ob=0, no sine). Free: D, mc; l0 FIXED=0. WLS vs Migkas 71 bins.
bubble_core has no explicit point-mass switch: CoE mass is exponential, M(<r)=M*Fx(r/LD)/Fx(1/LD).
Point mass = LD->0 limit; we fix LD=1e-4 (units of D) so Fx(r/LD)=1 to machine precision for all
r>~4e-3 D. D capped at 1e4 Gpc: beyond that the integrator is numerically unstable (spurious chi2 dips at D~1e15). sky() floor LD>0.01 bypassed by calling run directly."""
import json, importlib.util, numpy as np
from scipy.optimize import minimize
spec=importlib.util.spec_from_file_location("m","/workspace/cosmology-course/notes/catwise-sigma-null.py")
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
from bubble_core import run, shifted, IMIG, Fx
LD_PT=1e-4; CN=45.13662301218721
def pred(D,mc,l0):
    if not (m.RV*1.05<D<1e4 and 1e-9<mc<50): return None
    rr=run(D/m.RV,LD_PT,mc,0.0); af=rr["af"]
    if not np.all(np.isfinite(af)) or rr["vrmin"]<=0: return None
    return shifted(af/af.mean(),l0),rr
def obj(u):
    o=pred(10**u[0],10**u[1],0.0)
    return 800.0 if o is None else float(np.sum(m.W*(o[0][IMIG]-m.Y)**2))
starts=[(D,mc,l0) for D in (16,25,40,60,90) for mc in (1e-6,1e-5,1e-4,1e-3,3e-3) for l0 in (0,)]
res=[]
for s in starts:
    u0=np.array([np.log10(s[0]),np.log10(s[1])])
    if obj(u0)>=800: continue
    r=minimize(obj,u0,method="Nelder-Mead",options=dict(maxiter=500,xatol=1e-5,fatol=1e-6,adaptive=True))
    res.append((r.fun,r.x)); print(s,round(r.fun,4),10**r.x[0],10**r.x[1],0.0,flush=True)
res.sort(key=lambda t:t[0]); bc,bu=res[0]
r=minimize(obj,bu,method="Nelder-Mead",options=dict(maxiter=1000,xatol=1e-7,fatol=1e-8,adaptive=True))
if r.fun<bc: bc,bu=r.fun,r.x
D,mc,l0=10**bu[0],10**bu[1],0.0
full,rr=pred(D,mc,l0); p=full[IMIG]
# check point-mass limit insensitivity
ld_chk={ld:float(np.sum(m.W*(shifted((lambda a:a/a.mean())(run(D/m.RV,ld,mc,0.0)["af"]),l0)[IMIG]-m.Y)**2)) for ld in (1e-3,1e-4,1e-5)}
row=m.pack_row("pointmass_l0zero_D_mc",bc,2,["D","mc"],m.r2_u(m.Y,p),m.r2_w(m.Y,p,m.W),CN,
 extra=dict(D_Gpc=D,mc=mc,l0=l0,LD_fixed=LD_PT,Ob=0.0,shell_mass=0.0,sine=False,chi2_null=CN,
  vrmin=rr["vrmin"],Ht=rr["Ht"],age_Gyr=rr["age"],Fx_at_r0p01D=float(Fx(0.01/LD_PT)),chi2_vs_LD_check=ld_chk,
  n_starts=len(res),top5_chi2=[float(t[0]) for t in res[:5]],
  l0_fixed=True,note="chi2 flattens for D>~300 Gpc (far-field limit); D cap 1e4 Gpc. Point mass realised as LD->0 limit of bubble_core exponential CoE profile (LD=1e-4 D)."))
json.dump(row,open("/workspace/cosmology-course/notes/catwise-pointmass-l0zero.json","w"),indent=2,default=str)
print(json.dumps(row,indent=1,default=str))

"""Pure point mass at CoE (LD=1e-4, no shell, Ob=0, no sine), l0 fixed 0. Free: D, mc, R (obs distance), 0.1<=R<=14.3 Gpc.
NOTE: bubble_core.run depends on geometry only via DR=D/R (RKM enters only the c34 diagnostic), so D and R are degenerate:
chi2 is a function of (D/R, mc) only. Constraint D>1.05R, D<1e4 Gpc kept from the base script."""
import json, importlib.util, numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2 as c2, norm
spec=importlib.util.spec_from_file_location("m","/workspace/cosmology-course/notes/catwise-sigma-null.py")
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
from bubble_core import run, shifted, IMIG
LD_PT=1e-4; CN=45.13662301218721; RLO,RHI=0.1,14.3
def pred(D,mc,R):
    if not (RLO<=R<=RHI and R*1.05<D<1e4 and 1e-9<mc<50): return None
    rr=run(D/R,LD_PT,mc,0.0); af=rr["af"]
    if not np.all(np.isfinite(af)) or rr["vrmin"]<=0: return None
    return shifted(af/af.mean(),0.0),rr
def obj(u):
    o=pred(10**u[0],10**u[1],10**u[2])
    return 800.0 if o is None else float(np.sum(m.W*(o[0][IMIG]-m.Y)**2))
res=[]
for R in (1,3,6,10,14.3):
  for DR in (1.2,2,4,7,15):
    for mc in (1e-5,1e-4,1e-3,3e-3):
      u0=np.log10([DR*R,mc,R])
      if obj(u0)>=800: continue
      r=minimize(obj,u0,method="Nelder-Mead",options=dict(maxiter=800,xatol=1e-5,fatol=1e-6,adaptive=True))
      res.append((r.fun,r.x,R)); print(R,DR,mc,round(r.fun,4),10**r.x,flush=True)
res.sort(key=lambda t:t[0]); bc,bu,_=res[0]
r=minimize(obj,bu,method="Nelder-Mead",options=dict(maxiter=2000,xatol=1e-7,fatol=1e-8,adaptive=True))
if r.fun<bc: bc,bu=r.fun,r.x
D,mc,R=10**bu; full,rr=pred(D,mc,R); p=full[IMIG]
dchi=CN-bc; pv=c2.sf(dchi,3); sig=norm.isf(pv)
best_by_R={}
for t in res:
    k=str(t[2]); best_by_R[k]=min(best_by_R.get(k,1e9),float(t[0]))
row=dict(name="pointmass_l0zero_D_mc_Rfree",chi2=float(bc),n_free=3,free_params=["D","mc","R"],nu=int(m.N-3),
 R2_unweighted=float(m.r2_u(m.Y,p)),R2_weighted=float(m.r2_w(m.Y,p,m.W)),D_Gpc=float(D),mc=float(mc),R_Gpc=float(R),DR=float(D/R),
 chi2_null=CN,Delta_chi2=float(dchi),p_nested_nu3=float(pv),sigma_nested_nu3=float(sig),
 R_at_lower_bound=bool(R<=RLO*1.001),R_at_upper_bound=bool(R>=RHI*0.999),R_bounds=[RLO,RHI],
 best_chi2_by_start_R=best_by_R,top5_chi2=[float(t[0]) for t in res[:5]],n_starts=len(res),vrmin=rr["vrmin"],
 note="Model depends only on D/R and mc; R is fully degenerate with D (flat direction), so R value is arbitrary within bounds subject to D<1e4 Gpc.")
json.dump(row,open("/workspace/cosmology-course/notes/catwise-pointmass-l0zero-Rfree.json","w"),indent=2,default=str)
print(json.dumps(row,indent=1,default=str))

import json,sys
sys.argv=['x']
import importlib.util
spec=importlib.util.spec_from_file_location("m","/workspace/cosmology-course/notes/catwise-sigma-null.py")
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
import numpy as np
from scipy.optimize import minimize
prev=json.load(open("/workspace/cosmology-course/notes/catwise-sigma-null.json"))
CN=45.13662301218721; CF=8.4753
obj=lambda u: m.chi2_free(10**u[0],10**u[1],10**u[2],0.0)
starts=[(30,.15,.1),(40,.09,.04),(35,.12,.08),(45,.08,.03),(25,.2,.15),(50,.1,.05),(38,.0877,.0389),(28,.18,.2),(27.95,.217,.19),(20,.3,.3),(60,.06,.02)]
bc,bu=np.inf,None
for s in starts:
    r=minimize(obj,np.log10(s),method="Nelder-Mead",options=dict(maxiter=400,xatol=1e-4,fatol=1e-4,adaptive=True))
    print(s,r.fun,flush=True)
    if r.fun<bc: bc,bu=r.fun,r.x
r=minimize(obj,bu,method="Nelder-Mead",options=dict(maxiter=600,xatol=1e-5,fatol=1e-5,adaptive=True))
if r.fun<bc: bc,bu=r.fun,r.x
D,LD,mc=10**bu
full,rr=m.model_bins(D,LD,mc,0.0); p=full[m.IMIG]
row=m.pack_row("pure_exp_DL_mc_l0zero",bc,3,["D","L","mc"],m.r2_u(m.Y,p),m.r2_w(m.Y,p,m.W),CN,
 extra=dict(D_Gpc=D,LD=LD,L_Gpc=LD*D,mc=mc,l0_fixed=0.0,Ob=0.0,chi2_null=CN,chi2_l0free=CF,
 Delta_chi2_cost_vs_l0free=bc-CF,sigma_l0_vs_l0zero=m.sigma_from_delta(bc-CF,1)[0],vrmin=float(rr["vrmin"])))
json.dump(row,open("/workspace/cosmology-course/notes/catwise-exp-l0zero.json","w"),indent=2)
print(json.dumps(row,indent=1))

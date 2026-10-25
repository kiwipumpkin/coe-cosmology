exec(open('point_c25.py').read().split("# multistart")[0])
from scipy.optimize import minimize
import json
rng=np.random.default_rng(2);res={}
for uu in (2.5,5,10,30):
    lu=np.log10(uu);best=None
    for k in range(25):
        x0=LO+(HI-LO)*rng.random(7);x0=np.delete(x0,1)
        f=lambda y:-lnp(np.insert(y,1,lu)) if np.isfinite(lnp(np.insert(y,1,lu))) else 1e9
        r=minimize(f,x0,method='Nelder-Mead',options=dict(maxiter=6000,xatol=1e-6,fatol=1e-6))
        if best is None or r.fun<best.fun: best=r
    x=np.insert(best.x,1,lu);m=model(x);res[uu]=dict(x=x.tolist(),chi2=2*best.fun,R2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)))
    print(uu,round(2*best.fun,4),round(res[uu]['R2'],4),np.round(x,4).tolist())
json.dump(res,open('prof_c25.json','w'),indent=1)

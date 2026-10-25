import json,numpy as np
from scipy.optimize import differential_evolution, minimize
import joint_vgtc as J
B=J.B
def ev(x):
    if any(v<lo or v>hi for v,(lo,hi) in zip(x,B)): return None
    try: return J.evalp(x)
    except Exception: return None
def fchi(x):
    e=ev(x); return 1e6 if e is None else e['chi2']
def frms(x,T):
    e=ev(x); return 1e6 if e is None else e['rms']+1.0*max(0,e['chi2']-T)
def run(f,args=(),x0=None):
    best=None
    for sd in (1,2):
        r=differential_evolution(f,B,args=args,seed=sd,popsize=25,maxiter=500,tol=1e-10,polish=False,x0=x0)
        r2=minimize(f,r.x,args=args,method='Nelder-Mead',options=dict(maxiter=15000,xatol=1e-10,fatol=1e-13,adaptive=True))
        x,v=(r2.x,r2.fun) if r2.fun<r.fun else (r.x,r.fun)
        if best is None or v<best[1]: best=(x,v)
    return best
out={}
xc,cmin=run(fchi); e=ev(xc); out['min_chi2']=dict(x=list(xc),chi2=e['chi2'],cancel=e['cancel'],floor=e['rms']/1e-5,age=float(e['age']),Hmean=e['Hmean'],vo=e['vo'],M=e['M'],D=float(e['D']))
print('minchi2',out['min_chi2'],flush=True)
front=[];x0=xc
for T in [11.6,15,20,25,30,35,40,45]:
    if T<cmin-1e-9: front.append(dict(T=T,feasible=False)); continue
    x,v=run(frms,(T,),x0=x0); e=ev(x)
    d=dict(T=T,feasible=e['chi2']<=T*1.001,x=list(x),chi2=e['chi2'],cancel=e['cancel'],floor=e['rms']/1e-5,age=float(e['age']),Hmean=e['Hmean'],vo=e['vo'],M=e['M'],D=float(e['D']))
    front.append(d); print(T,{k:d[k] for k in d if k!='x'},flush=True)
out['front']=front
json.dump(out,open('pareto.json','w'),indent=1,default=float)

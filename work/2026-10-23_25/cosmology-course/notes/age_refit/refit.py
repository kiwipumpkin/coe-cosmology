import sys,json,numpy as np
sys.path.insert(0,'..'); import importlib.util
sp=importlib.util.spec_from_file_location('cc','../cmb-cancel.py'); cc=importlib.util.module_from_spec(sp); sp.loader.exec_module(cc)
from scipy.integrate import quad
from scipy.optimize import differential_evolution, minimize
HUB=977.8/70  # 1/H0 in Gyr at H0=70
def F(b1,mu,s):
    vo2=b1**2+2*mu*(1/s-1)
    if vo2<=0 or vo2>=1: return None
    vo=np.sqrt(vo2); Uo2=vo2/(1-vo2)
    f=lambda r: np.sqrt(r)/np.sqrt(Uo2*r+2*mu*(1-r/s))
    I=quad(f,0,s,limit=200)[0]
    return vo/s*I, vo
def obj(p,need):
    b1,q,mu,s,phi=p
    o=cc.obj(p,cc.G_FULL,2,B)
    a=F(b1,mu,s)
    if a is None: return 1.0
    return o+ 10*max(0,need-a[0])
B=[(1e-4,0.999),(1e-4,3.0),(0.0,2.0),(1.0,3.0),(0.0,360.0)]
base=np.array([0.5865,0.2304,0.09327,3.0,15.16])
out={}
def rep(x,name):
    b1,q,mu,s,phi=x; Fv,vo=F(b1,mu,s)
    S,Dm,Dg=cc.S_model(*x); R=cc.resid(cc.G_FULL,S)
    S0=cc.resid(cc.G_FULL,np.ones(360))
    d=dict(beta1=b1,q=q,mu=mu,s=s,phi=phi,v_o=vo,ageH=Fv,age_at_H70=Fv*HUB,H_needed_for_13=Fv*977.8/13,
           rms=cc.rms(R),floor_units=cc.rms(R)/1e-5,cancel_pct=100*(1-cc.rms(R)/cc.rms(S0)))
    out[name]=d; print(name,json.dumps(d)); return d
rep(base,'baseline')
# ceiling: max F over feasible box ignoring CMB
best=None
for sd in (1,2,3):
    r=differential_evolution(lambda p,need=0.9306: obj(p,need),B,seed=sd,popsize=30,maxiter=1500,tol=1e-10,polish=False,x0=base)
    r2=minimize(lambda p: obj(p,0.9306),r.x,method='Nelder-Mead',options=dict(xatol=1e-10,fatol=1e-14,maxiter=20000,adaptive=True))
    x=r2.x if r2.fun<r.fun else r.x
    f=min(r.fun,r2.fun); print('seed',sd,f,x)
    if best is None or f<best[1]: best=(x,f)
rep(best[0],'constrained_H70_age13')
json.dump(out,open('refit.json','w'),indent=1,default=float)

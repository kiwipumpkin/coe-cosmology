import sys,json,numpy as np; sys.path.insert(0,'/workspace/cosmology-course/notes/paper_repro/grfit'); import warnings; warnings.filterwarnings('ignore')
from core import *
def par(u): return 10**u[0],10**u[1],1+10**u[2],1+10**u[3],u[4]
def ev(u):
    if not (12<=u[4]<=18): return 1e9
    D,mcc,kin,kout,l0=par(u); af=ltb(D,mcc,kin,kout)
    if af is None or not np.all(np.isfinite(af)): return 1e9
    return chi(shifted(af,l0))
if __name__=='__main__':
    b=json.load(open('/workspace/cosmology-course/notes/paper_repro/grfit/fit_ltb_plain.json'))['rows'][4]
    u0=np.array([np.log10(b['D']),np.log10(b['mcc']),np.log10(b['kin']-1),np.log10(b['kout']-1),b['l0']])
    from scipy.optimize import minimize
    r=minimize(ev,u0,method='Nelder-Mead',options=dict(xatol=1e-7,fatol=1e-8,maxiter=4000,adaptive=True)); u0=r.x; c0=r.fun; print('best',par(u0),c0,flush=True)
    rng=np.random.default_rng(0); step=np.array([0.02,0.02,0.05,0.05,0.3]); u=u0.copy(); c=c0; S=[]; acc=0
    for i in range(40000):
        v=u+step*rng.standard_normal(5); cv=ev(v)
        if np.log(rng.random())< -(cv-c)/2/2.0: u,c=v,cv; acc+=1
        if c<c0: u0,c0=u.copy(),c
        S.append([*u,c])
    S=np.array(S); print('accept',acc/len(S),'best',par(u0),c0)
    np.save('chain.npy',S); json.dump(dict(u0=u0.tolist(),chi0=c0),open('best.json','w'))

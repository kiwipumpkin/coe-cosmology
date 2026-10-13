"""Multi-start fits for one configuration class. usage: python3 fit_run.py KINDS NSTARTS SEED -> raw/KINDS.json"""
import sys, os, json, time, numpy as np
from scipy.optimize import least_squares
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from ss_core import *
SC=1/(FLOOR*np.sqrt(360))
def fit_one(m,p0):
    lam=30.0; p=p0
    for it in range(4):
        f=lambda q: np.concatenate([(R-m.pattern(q))*SC, lam*m.violations(q)])
        r=least_squares(f,p,bounds=(m.lo,m.hi),method="trf",x_scale="jac",max_nfev=60*len(p),xtol=1e-8,ftol=1e-8)
        p=r.x
        if m.violations(p).max(initial=0)<=1e-3: break
        lam*=10
    return p
def main(kinds,ns,seed):
    m=Model(kinds); rng=np.random.default_rng(seed); out=[]; t0=time.time()
    for i in range(ns):
        p0=np.zeros(4*m.n).reshape(m.n,4)
        p0[:,0]=rng.uniform(0,0.6,m.n); p0[:,1]=rng.uniform(-np.pi,np.pi,m.n); p0[:,2]=np.arcsin(rng.uniform(-1,1,m.n)); p0[:,3]=rng.uniform(0,1,m.n)
        p=fit_one(m,p0.ravel()); X,M,Rv,Rc=m.decode(p); w=m.pattern(p)
        out.append(dict(left=float(rms(R-w)/FLOOR),viol=float(m.violations(p).max(initial=0)),X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),Rc=Rc.tolist(),
                        amps=(amps(R-w)/FLOOR).tolist()))
    os.makedirs("raw",exist_ok=True)
    json.dump(dict(kinds=kinds,nstarts=ns,seed=seed,secs=time.time()-t0,runs=out),open(f"raw/{kinds}.json","w"))
    L=sorted(r["left"] for r in out if r["viol"]<=1e-3)
    print(kinds,f"{time.time()-t0:.0f}s best={L[0] if L else None:.1f} med={np.median(L) if L else None:.1f} feas={len(L)}/{ns}",flush=True)
if __name__=="__main__": main(sys.argv[1],int(sys.argv[2]),int(sys.argv[3]))

import os, sys, json, numpy as np
os.environ["OMP_NUM_THREADS"]="1"
from multiprocessing import Pool
from mix_core import *
CASES={"inside_mirror":lambda:[G(0.5,0.5,13.8,0)],"both_nomirror":lambda:[G(0.5,0.5,13.8,1),G(0.6,14.8,30.0,-1)]}
_g={}
def job(a):
    case,X,F,lam=a
    if case not in _g: _g[case]=CASES[case]()
    gs=_g[case]; ds,L=solve_box(gs,X,F,lam)
    d=np.concatenate(ds); vol=np.concatenate([np.full(len(g.P),g.H**3) for g in gs]); m=np.concatenate([g.m for g in gs])
    np.savez(f"sol_{case}_X{X}_F{F}_lam{lam}_it{os.environ.get("MAXIT","400")}.npz",*ds)
    lo,hi=bounds(X,F)
    return dict(case=case,X=X,F=F,lam=lam,maxit=os.environ.get("MAXIT","400"),left=L,M_abs=float(np.abs(m*d).sum()),peak=float(np.abs(d).max()),
        vol_active=float(vol[np.abs(d)>0.005].sum()),vol_at_bound=float(vol[(d>=0.98*hi)|(d<=0.98*lo)].sum()))
if __name__=="__main__":
    tasks=[(c,X,float(sys.argv[2]) if len(sys.argv)>2 else 1.0,float(sys.argv[3]) if len(sys.argv)>3 else 0.0) for c,Xs in json.loads(sys.argv[1]).items() for X in Xs]
    out=[]
    with Pool(4) as p:
        for r in p.imap_unordered(job,tasks):
            print(json.dumps(r),flush=True); out.append(r)

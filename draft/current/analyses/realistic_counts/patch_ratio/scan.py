import os, sys, json, numpy as np
from multiprocessing import Pool
from ps_core import *
_g={}
def job(a):
    c,B=a
    if c not in _g: _g[c]=CONFIGS[c]()
    ds,L=solve_sym(_g[c],B); np.savez(f"sol_{c}_B{B}.npz",*ds)
    d=np.concatenate(ds); vol=np.concatenate([np.full(len(g.P),g.H**3) for g in _g[c]])
    return dict(cfg=c,B=B,left=L,vol_active=float(vol[np.abs(d)>0.1*B].sum()),vol_at_bound=float(vol[np.abs(d)>=0.98*B].sum()))
if __name__=="__main__":
    tasks=[(c,B) for c,Bs in json.loads(sys.argv[1]).items() for B in Bs]
    with Pool(int(os.environ.get("NP","3"))) as p:
        for r in p.imap_unordered(job,tasks): print(json.dumps(r),flush=True)

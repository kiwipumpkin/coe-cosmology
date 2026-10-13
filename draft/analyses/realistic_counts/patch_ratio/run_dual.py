import os, sys, json, time, numpy as np
from multiprocessing import Pool
from ps_core import *
_g={}
def job(a):
    c,B=a; t=time.time()
    if c not in _g: _g[c]=CONFIGS[c]()
    ds,lp,ld=solve_dual(_g[c],B)
    return dict(cfg=c,B=B,left_primal=lp,left_dual=ld,sec=round(time.time()-t))
if __name__=="__main__":
    tasks=[(c,B) for c,Bs in json.loads(sys.argv[1]).items() for B in Bs]
    with Pool(int(os.environ.get("NP","2"))) as p:
        for r in p.imap_unordered(job,tasks): print(json.dumps(r),flush=True)

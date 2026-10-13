"""Cold random starts at a given composition (parallel). usage: python3 cold.py KINDS NSTARTS -> cold_KINDS.json"""
import os, sys, json, time, numpy as np
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from extend import job
def p0(n,seed):
    rng=np.random.default_rng(seed); p=np.zeros((n,4))
    p[:,0]=rng.uniform(0,0.6,n); p[:,1]=rng.uniform(-np.pi,np.pi,n); p[:,2]=np.arcsin(rng.uniform(-1,1,n)); p[:,3]=rng.uniform(0,1,n); return p.ravel()
if __name__=="__main__":
    k=sys.argv[1]; ns=int(sys.argv[2]); t0=time.time()
    with Pool(8) as pool: out=pool.map(job,[(k,p0(len(k),5000+i),i) for i in range(ns)])
    json.dump(out,open(f"cold_{k}.json","w"))
    L=sorted(r["left"] for r in out if r["viol"]<=1e-3)
    print(k,f"{time.time()-t0:.0f}s",np.round(L,2).tolist(),flush=True)

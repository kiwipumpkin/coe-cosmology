"""Short replacement for the slow bisection: dual at a few B in the known bracket, then a polished primal at the first B reaching 1x."""
import sys, json, numpy as np
from ps_core import *
c=sys.argv[1]; Bs=[float(x) for x in sys.argv[2].split(",")]
gs=CONFIGS[c](); A=np.concatenate([g.A for g in gs],1); m=np.concatenate([g.m for g in gs])/1e50; Ap=np.r_[A,1e3*m[None]]; y=np.r_[yR,0.0]
res=dict(cfg=c,dual={}); u=None; below=None
for B in Bs:
    uu,ld=dual_u(Ap,y,B,u); res["dual"][B]=ld; print(c,B,ld,flush=True)
    if ld>1: below=B; u=uu; continue
    d,L=min((polish(A,Ap,y,uu,B,nf) for nf in (100,400,1500)),key=lambda t:t[1]); print(" primal",L,flush=True)
    if L<=1:
        res.update(B_below=below,B_star=B,B_1x_primal=B,left_primal_1x=L); out=[];i=0
        for g in gs: out.append(d[i:i+len(g.P)]); i+=len(g.P)
        np.savez(f"thr_{c}.npz",*out); break
json.dump(res,open(f"thr_{c}.json","w"),indent=1); print(json.dumps(res))

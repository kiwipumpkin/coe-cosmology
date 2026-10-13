import sys, json, numpy as np
from ps_core import *
c=sys.argv[1]; lo,hi=float(sys.argv[2]),float(sys.argv[3])
gs=CONFIGS[c](); r,ds=threshold(gs,lo,hi); r["cfg"]=c
if r["left_primal_at_Bstar"]>1:
    A=np.concatenate([g.A for g in gs],1); m=np.concatenate([g.m for g in gs])/1e50; Ap=np.r_[A,1e3*m[None]]; y=np.r_[yR,0.0]
    for fac in (1.01,1.02,1.04):
        B=r["B_star"]*fac; u,ld=dual_u(Ap,y,B); d,L=min((polish(A,Ap,y,u,B,nf) for nf in (100,400,1500)),key=lambda t:t[1])
        r["B_1x_primal"]=B; r["left_primal_1x"]=L
        if L<=1: break
    out=[];i=0
    for g in gs: out.append(d[i:i+len(g.P)]); i+=len(g.P)
    ds=out; print(json.dumps(r),flush=True)
np.savez(f"thr_{c}.npz",*ds); json.dump(r,open(f"thr_{c}.json","w"))

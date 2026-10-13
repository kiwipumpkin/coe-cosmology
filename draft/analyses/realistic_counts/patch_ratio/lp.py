"""Minimal box size B for (a) full cancellation A d = yR, and (b) a guaranteed <=1x leftover (each residual |r_i|<=rtol/sqrt(n), inside the
L2 ball) — LP: min s s.t. -s<=d<=s, zero net mass.  Gives feasible primal solutions (upper bounds on the threshold)."""
import os, sys, json, time, numpy as np
from scipy.optimize import linprog
from scipy import sparse
from ps_core import *
def minB(grids,rtol):
    A=np.concatenate([g.A for g in grids],1); m=np.concatenate([g.m for g in grids])/1e50; n,N=A.shape
    # vars: d (N), s (1). minimize s.
    c=np.r_[np.zeros(N),1.0]
    I=sparse.identity(N,format="csr"); one=sparse.csr_matrix(np.ones((N,1)))
    Aub=[sparse.hstack([I,-one]),sparse.hstack([-I,-one])]; bub=[np.zeros(N),np.zeros(N)]
    tol=np.sqrt(max(rtol**2-TAIL**2,0)/n)
    if tol>0:
        Ad=sparse.csr_matrix(A); Aub+= [sparse.hstack([Ad,sparse.csr_matrix((n,1))]),sparse.hstack([-Ad,sparse.csr_matrix((n,1))])]
        bub+=[yR+tol,-yR+tol]; Aeq=sparse.csr_matrix(np.r_[m,0.0][None]); beq=[0.0]
    else:
        Aeq=sparse.csr_matrix(np.c_[np.r_[A,m[None]],np.zeros((n+1,1))]); beq=np.r_[yR,0.0]
    r=linprog(c,A_ub=sparse.vstack(Aub).tocsr(),b_ub=np.concatenate(bub),A_eq=Aeq,b_eq=beq,bounds=[(None,None)]*N+[(0,None)],method="highs")
    d=r.x[:N]; left=float(np.sqrt(((yR-A@d)**2).sum()+TAIL**2)); out=[];i=0
    for g in grids: out.append(d[i:i+len(g.P)]); i+=len(g.P)
    return float(r.x[-1]),left,out,r.status
if __name__=="__main__":
    res={}
    for c in sys.argv[1:]:
        t=time.time(); gs=CONFIGS[c](); res[c]={}
        for name,rt in [("exact",0.0),("1x",1.0)]:
            B,L,ds,st=minB(gs,rt); np.savez(f"lp_{c}_{name}.npz",*ds)
            res[c][name]=dict(B=B,left=L,status=st); print(c,name,B,L,st,round(time.time()-t),"s",flush=True)
        json.dump(res,open(f"lp_{'_'.join(sys.argv[1:])}.json","w"),indent=1)

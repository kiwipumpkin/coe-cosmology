"""Mixed cluster+void composition. Bookkeeping from ../../local_web/volumes_counts.md (foam PACK):
per web cell s^3 (s=0.05 Gpc): 1 void (deficit D_v = 0.546 rho s^3, volume 0.78 s^3) and 6.77 clusters, each carrying the positive
mass E_c = 0.546 rho s^3 / 6.77 (cluster 4.5e44 kg + its share of filaments/walls). So the cosmic mean ratio r0 = 6.77 is exactly net zero.
In a cell with cluster count x(1+a) and void count x(1+b):   delta = 0.546 (a - b),   ratio q = r/r0 = (1+a)/(1+b).
Ratio constraint |log q| <= log(1+X) and void filling 0.78(1+b) <= f_max give the exact delta bounds
     -0.546 F X/(1+X) <= delta <= 0.546 F X,     F = f_max/0.78   (main F=1: voids may not fill more than the mean 78%)."""
import os, sys, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.dirname(HERE)); _c=os.getcwd()
from field_multi import G, yR, TAIL, MM, rho_gpc, GPC
os.chdir(_c)
from scipy.optimize import lsq_linear
from scipy.sparse.linalg import LinearOperator
S_GPC=0.05; RHO=2.9e-27; DV_FRAC=0.546; NC_PER_CELL=6.77; FV0=0.78
n_v0=1/S_GPC**3; n_c0=NC_PER_CELL*n_v0                     # per Gpc^3: 8000 voids, 54160 clusters
def bounds(X,F=1.0): return -DV_FRAC*F*X/(1+X), DV_FRAC*F*X
def solve_box(grids,X,F=1.0,lam=0.0,kappa=1e3):
    """min ||yR - A d||^2 + kappa^2 (mass)^2 + lam * sum(vol d^2), lo<=d<=hi.  Returns list of per-grid deltas and leftover."""
    A=np.concatenate([g.A for g in grids],1); m=np.concatenate([g.m for g in grids])/1e50; vol=np.concatenate([np.full(len(g.P),g.H**3) for g in grids])
    nd=A.shape[0]; N=A.shape[1]; sl=np.sqrt(lam*vol)
    def mv(x): return np.r_[A@x,kappa*(m@x),sl*x]
    def rmv(y): return A.T@y[:nd]+kappa*m*y[nd]+sl*y[nd+1:]
    op=LinearOperator((nd+1+N,N),matvec=mv,rmatvec=rmv,dtype=float)
    lo,hi=bounds(X,F); b=np.r_[yR,0.0,np.zeros(N)]
    r=lsq_linear(op,b,bounds=(np.full(N,lo),np.full(N,hi)),method="trf",lsmr_tol="auto",max_iter=int(os.environ.get("MAXIT","400")),tol=1e-10)
    d=r.x; left=float(np.sqrt(((yR-A@d)**2).sum()+TAIL**2)); out=[];i=0
    for g in grids: out.append(d[i:i+len(g.P)]); i+=len(g.P)
    return out,left
def composition(d):
    """conversion rule (total object count fixed: clusters <-> voids):  a*n_c0 = -b*n_v0  -> delta = 0.546*a*(1+6.77) ; q=(1+a)/(1+b)"""
    a=d/(DV_FRAC*(1+NC_PER_CELL)); b=-NC_PER_CELL*a; q=(1+a)/(1+b)
    return a,b,q

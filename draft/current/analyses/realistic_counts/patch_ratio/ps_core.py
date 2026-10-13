"""Patch-ratio model. A domain (grid cell, size H) = many patches of scale r_p; over-patches have d_plus, under-patches d_minus.
Over-patch volume fraction f; cosmic mean f0 = -d_minus/(d_plus-d_minus) gives net zero.  Domain contrast d = (f-f0)*(d_plus-d_minus).
So |f-f0|<=X  <=>  |d| <= B = X*Delta  : a symmetric box problem in d (solved here), mapped to f per patch choice afterwards."""
import os, sys, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.dirname(HERE)); _c=os.getcwd()
from field_multi import G, yR, TAIL
os.chdir(_c)
from scipy.optimize import lsq_linear
from scipy.sparse.linalg import LinearOperator
CONFIGS={"in05":lambda:[G(0.5,0.5,13.8,0)],"in10":lambda:[G(1.0,0.5,13.8,0)],"in20":lambda:[G(2.0,0.5,13.8,0)],
         "both05":lambda:[G(0.5,0.5,13.8,1),G(0.6,14.8,30.0,-1)],"both10":lambda:[G(1.0,0.5,13.8,1),G(1.0,14.8,30.0,-1)]}
# patch choices: r_p = fit radius scale (Mpc); values from ../../density_vs_scale/fits_v2.json (over_noDef, under_noDisp power laws)
PATCH={"extreme 30 Mpc":(30,2.10,-0.55),"extreme 50 Mpc":(50,0.83,-0.39),"extreme 100 Mpc":(100,0.24,-0.25),
       "typical ±0.3, 50 Mpc":(50,0.30,-0.30),"typical ±0.3, 30 Mpc":(30,0.30,-0.30),"typical ±0.3, 100 Mpc":(100,0.30,-0.30)}
def pinfo(name):
    rp,dp,dm=PATCH[name]; D=dp-dm; f0=-dm/D; return dict(rp=rp,dp=dp,dm=dm,Delta=D,f0=f0,Vp=4/3*np.pi*(rp/1000)**3)
def solve_sym(grids,B,kappa=1e3):
    A=np.concatenate([g.A for g in grids],1); m=np.concatenate([g.m for g in grids])/1e50; nd=A.shape[0]; N=A.shape[1]
    op=LinearOperator((nd+1,N),matvec=lambda x:np.r_[A@x,kappa*(m@x)],rmatvec=lambda y:A.T@y[:nd]+kappa*m*y[nd],dtype=float)
    r=lsq_linear(op,np.r_[yR,0.0],bounds=(np.full(N,-B),np.full(N,B)),method="trf",lsmr_tol="auto",max_iter=int(os.environ.get("MAXIT","1500")),tol=1e-10)
    d=r.x; left=float(np.sqrt(((yR-A@d)**2).sum()+TAIL**2)); out=[];i=0
    for g in grids: out.append(d[i:i+len(g.P)]); i+=len(g.P)
    return out,left
from scipy.optimize import minimize
def solve_dual(grids,B,kappa=1e3,eps_list=(3.0,1.0,0.3,1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5)):
    """Exact-ish box LSQ via the 41-dim dual: max_u y.u - |u|^2/2 - B*sum_j |a_j.u| (smoothed |z|~sqrt(z^2+e^2), e relative to
    typical |a_j.u|, continuation).  Primal d_j = B*v_j/sqrt(v_j^2+e^2) is always inside the box; residual = u at optimum.
    Returns (ds, leftover_primal(upper bound), leftover_dual(lower bound))."""
    A=np.concatenate([g.A for g in grids],1); m=np.concatenate([g.m for g in grids])/1e50
    Ap=np.r_[A,kappa*m[None]]; y=np.r_[yR,0.0]; u=y.copy()
    best=(np.inf,None)
    for e_rel in eps_list:
        sc=np.median(np.abs(Ap.T@u))+1e-30; e=e_rel*sc
        def f(u):
            v=Ap.T@u; s=np.sqrt(v*v+e*e); val=y@u-0.5*u@u-B*s.sum()
            g=y-u-B*(Ap@(v/s)); return -val,-g
        u=minimize(f,u,jac=True,method="L-BFGS-B",options=dict(maxiter=5000,gtol=1e-12,ftol=1e-15)).x
        v=Ap.T@u; dd=B*v/np.sqrt(v*v+e*e); res=yR-A@dd; ob=(res**2).sum()+(kappa*(m@dd))**2
        if ob<best[0]: best=(ob,dd)
    d=best[1]; res=yR-A@d; lp=float(np.sqrt((res**2).sum()+TAIL**2)); v=Ap.T@u
    dual=y@u-0.5*u@u-B*np.abs(v).sum(); ld=float(np.sqrt(max(2*dual,0)+TAIL**2))  # lower bound on sqrt(|res|^2+mass term)
    out=[];i=0
    for g in grids: out.append(d[i:i+len(g.P)]); i+=len(g.P)
    return out,lp,ld
def dual_u(Ap,y,B,u0=None,eps_list=(3.0,1.0,0.3,1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4,3e-5,1e-5)):
    u=y.copy() if u0 is None else u0.copy()
    for e_rel in eps_list:
        e=e_rel*(np.median(np.abs(Ap.T@u))+1e-30)
        def f(u):
            v=Ap.T@u; s=np.sqrt(v*v+e*e); return -(y@u-0.5*u@u-B*s.sum()),-(y-u-B*(Ap@(v/s)))
        u=minimize(f,u,jac=True,method="L-BFGS-B",options=dict(maxiter=5000,gtol=1e-12,ftol=1e-15)).x
    v=Ap.T@u; dual=y@u-0.5*u@u-B*np.abs(v).sum(); return u,float(np.sqrt(max(2*dual,0)+TAIL**2))
def polish(A,Ap,y,u,B,nfree=400):
    """complementary slackness: cells with |a_j.u| clearly >0 sit at +-B (sign of a_j.u); re-solve the few smallest-|v| cells freely in the box."""
    v=Ap.T@u; d=B*np.sign(v); fr=np.argsort(np.abs(v))[:nfree]; d[fr]=0
    r=lsq_linear(Ap[:,fr],y-Ap@d,bounds=(-B,B),method="bvls",tol=1e-12); d[fr]=r.x
    return d,float(np.sqrt(((yR-A@d)**2).sum()+TAIL**2))
def threshold(grids,Blo,Bhi,kappa=1e3,it=12):
    """bisection on the dual-certified minimal leftover = 1x; returns B*, polished primal at B* (and its leftover)."""
    A=np.concatenate([g.A for g in grids],1); m=np.concatenate([g.m for g in grids])/1e50; Ap=np.r_[A,kappa*m[None]]; y=np.r_[yR,0.0]
    u=None
    for _ in range(it):
        Bm=0.5*(Blo+Bhi); uu,ld=dual_u(Ap,y,Bm,u)
        if ld>1.0: Blo=Bm; u=uu
        else: Bhi=Bm
    u,ld=dual_u(Ap,y,Bhi,u); best=None
    for nf in (100,400,1500):
        d,L=polish(A,Ap,y,u,Bhi,nf)
        if best is None or L<best[1]: best=(d,L)
    d,L=best; out=[];i=0
    for g in grids: out.append(d[i:i+len(g.P)]); i+=len(g.P)
    return dict(B_star=Bhi,B_below=Blo,left_dual_at_Bstar=ld,left_primal_at_Bstar=L),out

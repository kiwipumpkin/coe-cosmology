"""Shallow, observation-based voids (Clampitt & Jain 2015, arXiv:1404.1834) for the sequential runs.
Profile (spherically symmetric, so the tidal pull at any point is exact from the enclosed mass M_enc(r)):
  core  r<Rv : delta(r) = -0.5 + 0.25 (r/Rv)^2      -> centre -0.5, mean inside Rv = -0.35, edge -0.25
  shell Rv..3Rv: compensating mass, M_enc(r) = -M_def * (1 - S(x)),  x=(r-Rv)/(2Rv), S = 3x^2-2x^3 (smoothstep)
               -> still half under-compensated at 2Rv, fully compensated at 3Rv, net mass exactly zero per void
  r>3Rv : M_enc = 0  -> the void exerts NO pull outside 3Rv (Newton's shell theorem)
M_def = 0.35 * rho(r_CoE) * (4/3) pi Rv^3.
Tidal pattern: same formula as web_core._pattern (point mass M -> M_enc(r) evaluated at each ring point and at the observer)."""
import os, sys, numpy as np, numba as nb
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.dirname(HERE))
_c=os.getcwd()
from seq_core import *
from web_core import TH
os.chdir(_c)
DELTA_MEAN=-0.35

def def_mass(x,Rv):   # kg, positive number = mass missing from the core
    return 0.35*rho_gpc(np.atleast_2d(x))*vol(np.asarray(Rv))*GPC**3

@nb.njit(cache=True)
def menc_frac(r,Rv):
    """M_enc(r)/M_def (negative inside 3Rv, 0 outside)."""
    if r>=3*Rv: return 0.0
    if r<=Rv:
        u=r/Rv   # integral of delta(r) 4 pi r^2 dr / (0.35 * 4/3 pi Rv^3), delta = -0.5+0.25u^2
        return (-0.5*u**3+0.25*0.6*u**5)/0.35
    x=(r-Rv)/(2*Rv); return -(1.0-(3*x*x-2*x*x*x))

@nb.njit(parallel=True,cache=True)
def _pat(X,Y,Z,K,Rv,cth,sth):
    n=X.shape[0]; nt=cth.shape[0]; out=np.zeros(nt)
    for t in nb.prange(nt):
        x=cth[t]; y=sth[t]; acc=0.0
        for i in range(n):
            dx=X[i]-x; dy=Y[i]-y; dz=Z[i]
            r=np.sqrt(dx*dx+dy*dy+dz*dz); D=np.sqrt(X[i]*X[i]+Y[i]*Y[i]+Z[i]*Z[i])
            fr=menc_frac(r,Rv[i]); fD=menc_frac(D,Rv[i])
            if fr!=0.0: acc+=K[i]*fr*(dx*x+dy*y)/(r*r*r)
            if fD!=0.0: acc-=K[i]*fD*(X[i]*x+Y[i]*y)/(D*D*D)
        out[t]=acc
    return out-out.mean()

def shallow_pattern(Xg,Rvg,Mdef=None):
    """Xg (n,3) Gpc, Rvg (n,) Gpc. Returns Delta T/T ring pattern (same units as R)."""
    Xg=np.atleast_2d(Xg); Rvg=np.broadcast_to(np.asarray(Rvg,float),(len(Xg),)).copy()
    if Mdef is None: Mdef=def_mass(Xg,Rvg)
    P=Xg/RV
    return _pat(P[:,0].copy(),P[:,1].copy(),P[:,2].copy(),np.ascontiguousarray(Mdef/KG),Rvg/RV,np.cos(TH),np.sin(TH))

@nb.njit(parallel=True,cache=True)
def _bank(X,Y,Z,K,Rv,cth,sth):
    n=X.shape[0]; nt=cth.shape[0]; out=np.zeros((n,nt))
    for i in nb.prange(n):
        D=np.sqrt(X[i]*X[i]+Y[i]*Y[i]+Z[i]*Z[i]); fD=menc_frac(D,Rv[i]); m=0.0
        for t in range(nt):
            x=cth[t]; y=sth[t]; dx=X[i]-x; dy=Y[i]-y; dz=Z[i]
            r=np.sqrt(dx*dx+dy*dy+dz*dz); fr=menc_frac(r,Rv[i]); a=0.0
            if fr!=0.0: a+=K[i]*fr*(dx*x+dy*y)/(r*r*r)
            if fD!=0.0: a-=K[i]*fD*(X[i]*x+Y[i]*y)/(D*D*D)
            out[i,t]=a; m+=a
        m/=nt
        for t in range(nt): out[i,t]-=m
    return out

def bank_patterns(Xg,Rv):
    """(n,360) pattern of each void separately (one batched call)."""
    Xg=np.atleast_2d(Xg); Rvg=np.broadcast_to(np.asarray(Rv,float),(len(Xg),)).copy(); M=def_mass(Xg,Rvg); P=Xg/RV
    return _bank(P[:,0].copy(),P[:,1].copy(),P[:,2].copy(),np.ascontiguousarray(M/KG),Rvg/RV,np.cos(TH),np.sin(TH))

if __name__=="__main__":
    # checks: (1) far from everything -> exactly zero; (2) point-mass limit: a deep 'void' with fr=-1 everywhere equals pattern() of -M_def
    x=np.array([[0,35.0,0]]); print("far void pattern rms/floor:",rms(shallow_pattern(x,1.5))/FLOOR)
    x=np.array([[0,16.5,0]]); print("near-ring void rms/floor:",rms(shallow_pattern(x,1.5))/FLOOR, " deep point-mass equiv:",rms(pattern(x,-def_mass(x,1.5)))/FLOOR)
    r=np.linspace(0,3.2,9); print("M_enc/M_def at r/Rv",[round(menc_frac(q,1.0),3) for q in r])

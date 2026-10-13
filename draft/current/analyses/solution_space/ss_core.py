"""Solution-space mapping: deliberately fit N concentrated masses and voids to cancel the CMB residual R.
Read-only reuse of ../local_web/web_core.py (R, FLOOR, KG, geometry, exp-lock density, field3 via lattice_core) and the
smooth net-mass compensation pattern of the four-clump baseline (../shell_ensemble/with_clumps/baseline_patterns.npz).

Object types (all positions free in 14.3 Gpc + core radius .. 43 Gpc from us, any direction, core outside the CoE interior):
  C  compact cluster : point mass M in [1e49, 2e52] kg (log-free); core = delta=+100 sphere (as the greedy near-us fit)
  V  deep void       : uncompensated, delta=-0.70, radius Rv in [0.5, 4] Gpc -> mass -0.7 rho V (point-mass pull outside the core)
  S  shallow void    : Clampitt & Jain profile (mean -0.35 inside Rv, compensation spread to 3 Rv, net zero), Rv in [0.5, 4] Gpc
Net mass: any non-zero net (C and V) is balanced by the same smooth removal/addition over the surround region that the
four-clump baseline used (pattern Fcomp_unit), so every fit has net mass zero, exactly as the 88.8x four-clump reference.
Cores never overlap (|xi-xj| >= Ri+Rj+0.02 Gpc)."""
import os, sys, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); CS=os.path.dirname(HERE)
sys.path.insert(0,os.path.join(CS,"local_web")); _c=os.getcwd()
from web_core import R,FLOOR,KG,RV,rms,C_G,A_IN_G,GPC,RHO,RMIN,rho_gpc,TH
import lattice_core as lc
os.chdir(_c)
FC=np.load(os.path.join(CS,"shell_ensemble","with_clumps","baseline_patterns.npz"))["Fcomp_unit"]
RMAX=43.0; GAP=0.02; DELTA_C=100.0
MLO,MHI=1e49,2e52; RLO,RHI=0.5,4.0
CTH,STH=np.cos(TH),np.sin(TH)
AXIS=15.16   # residual is a pure cosine series about l=15.16 / 195.16 (CoE direction): mirror symmetry of the problem

def vol(r): return 4/3*np.pi*r**3

def unit_field(X):
    """(n,360) pattern per 1 kg of point mass at X (Gpc), identical to web_core.pattern / field3."""
    P=np.atleast_2d(X)/RV; D=np.linalg.norm(P,axis=1)
    return lc.field3(np.arctan2(P[:,1],P[:,0]),np.arcsin(np.clip(P[:,2]/D,-1,1)),D)/KG

def menc_frac(r,Rv):
    u=r/Rv; x=(r-Rv)/(2*Rv)
    core=(-0.5*u**3+0.15*u**5)/0.35
    shell=-(1.0-(3*x*x-2*x**3))
    return np.where(r>=3*Rv,0.0,np.where(r<=Rv,core,shell))

def shallow_field(X,Rv):
    """(n,360) pattern of shallow compensated voids (same formula as ../sequential/shallow/shallow_core._pat)."""
    X=np.atleast_2d(X)/RV; Rv=np.atleast_1d(Rv)/RV
    Mdef=0.35*rho_gpc(np.atleast_2d(X*RV))*vol(Rv*RV)*GPC**3/KG
    x=CTH[None]; y=STH[None]; Xx,Y,Z=X[:,0:1],X[:,1:2],X[:,2:3]
    dx,dy,dz=Xx-x,Y-y,Z; r=np.sqrt(dx*dx+dy*dy+dz*dz); D=np.linalg.norm(X,axis=1)[:,None]
    fr=menc_frac(r,Rv[:,None]); fD=menc_frac(D,Rv[:,None])
    a=Mdef[:,None]*(fr*(dx*x+dy*y)/r**3-fD*(Xx*x+Y*y)/D**3)
    return a-a.mean(1,keepdims=True)

class Model:
    """kinds: string like 'CCVS'. Params per object: u (radial fraction), lon (rad), lat (rad), s (size fraction in [0,1])."""
    def __init__(self,kinds):
        self.k=np.array(list(kinds)); self.n=len(kinds)
        self.lo=np.tile([0,-4*np.pi,-np.pi/2,0],self.n); self.hi=np.tile([1,4*np.pi,np.pi/2,1],self.n)
    def decode(self,p):
        p=p.reshape(self.n,4); u,lon,lat,s=p.T
        M=np.zeros(self.n); Rv=np.zeros(self.n); Rc=np.zeros(self.n)
        isC=self.k=="C"
        M[isC]=10**(np.log10(MLO)+s[isC]*np.log10(MHI/MLO))
        Rv[~isC]=RLO+s[~isC]*(RHI-RLO)
        # core radii (cluster depends on local density -> use mean density first, tiny effect; refined below)
        dirs=np.c_[np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]
        Rc[~isC]=Rv[~isC]
        Rc[isC]=(3*M[isC]/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3)
        d=RMIN+Rc+u*(RMAX-RMIN-Rc)
        X=dirs*d[:,None]
        isV=self.k=="V"
        M[isV]=-0.70*rho_gpc(X[isV])*vol(Rv[isV])*GPC**3
        return X,M,Rv,Rc
    def pattern(self,p):
        X,M,Rv,Rc=self.decode(p); w=np.zeros(360)
        pm=self.k!="S"
        if pm.any():
            w+=M[pm]@unit_field(X[pm]); w-=(M[pm].sum()/KG)*FC        # net-mass compensation (smooth, as 4-clump baseline)
        if (~pm).any(): w+=shallow_field(X[~pm],Rv[~pm]).sum(0)
        return w
    def violations(self,p):
        X,M,Rv,Rc=self.decode(p)
        v=[np.maximum(0,A_IN_G+Rc-np.linalg.norm(X-C_G,axis=1))]
        if self.n>1:
            i,j=np.triu_indices(self.n,1); d=np.linalg.norm(X[i]-X[j],axis=1)
            v.append(np.maximum(0,Rc[i]+Rc[j]+GAP-d))
        return np.concatenate(v)
    def leftover(self,p): return rms(R-self.pattern(p))/FLOOR

def amps(v,M=12):
    m=np.arange(1,M+1)[:,None]
    return np.hypot(2*np.mean(v[None]*np.cos(m*TH),1),2*np.mean(v[None]*np.sin(m*TH),1))

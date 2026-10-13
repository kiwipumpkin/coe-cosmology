"""Sequential, strictly non-overlapping void/cluster placement (v2).
Read-only reuse of ../local_web/web_core.py (tidal pattern, residual R, R4, floor, geometry, exp-lock density).

Objects are spheres. Hard exclusion: |x_i - x_j| >= R_i + R_j + GAP for EVERY pair (no exceptions).
Void:    sphere radius Rv, density contrast delta=-0.70  -> mass deficit 0.70*rho*V  (negative mass)
Cluster: sphere at delta=+1 (same convention as the four clumps): R_c = (3 M / (4 pi rho))^(1/3)
Net mass is kept at zero by giving clusters exactly the mass the voids removed.
Note: hard spheres can't fill 78% of space (random sequential packing jams at ~38%, lattice ~52%); we report the void volume fraction reached.
"""
import os, sys, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(os.path.dirname(HERE),"local_web"))
_cwd=os.getcwd()
from web_core import (R,R4,BASE0,BASE4,FLOOR,KG,RV,rms,amps,pattern,C_G,A_IN_G,L_G,D_G,GPC,RHO,RMIN,rho_gpc,PACK)
os.chdir(_cwd)
DELTA_V=PACK["delta_void"]          # -0.70
S_MAX_G=3*RV                         # 42.9 Gpc: outer edge of the surround region used in all earlier runs
LOCAL_RMAX=30.0                      # near/tidal region 14.3-30 Gpc
GAP=0.02                             # Gpc, extra clearance between every pair of spheres

def vol(r): return 4/3*np.pi*r**3
def void_mass(x,Rv): return DELTA_V*rho_gpc(np.asarray(x))*vol(Rv)*GPC**3            # kg (negative)
def cluster_radius(M,x): return (3*abs(M)/(4*np.pi*rho_gpc(np.asarray(x))*GPC**3))**(1/3)   # delta=+1 sphere, Gpc

def region_ok(x,r,rmin=RMIN,rmax=S_MAX_G):
    """whole sphere outside the 14.3 Gpc vision sphere and outside the CoE interior (<43.9 Gpc from CoE); centre within rmax."""
    x=np.atleast_2d(x); r=np.atleast_1d(r); s=np.linalg.norm(x,axis=1)
    return (s-r>=max(rmin,RMIN)) & (s<=rmax) & (np.linalg.norm(x-C_G,axis=1)-r>=A_IN_G)

def ldiff(w,target=None):
    return rms((R if target is None else target)-w)/FLOOR

def leftover(w,target=None): return ldiff(w,target)

class Packer:
    """Grid hash of placed spheres; exact pairwise exclusion test."""
    def __init__(self,h=2.0):
        self.h=h; self.grid={}; self.X=[]; self.r=[]; self.M=[]; self.kind=[]; self.rmax=0.0
    def _key(self,x): return tuple(np.floor(np.asarray(x)/self.h).astype(int))
    def free(self,x,r,gap=GAP):
        reach=r+self.rmax+gap; k=self._key(x); n=int(np.ceil(reach/self.h))
        for i in range(k[0]-n,k[0]+n+1):
            for j in range(k[1]-n,k[1]+n+1):
                for l in range(k[2]-n,k[2]+n+1):
                    for idx in self.grid.get((i,j,l),()):
                        d=self.X[idx]-x
                        if d@d < (r+self.r[idx]+gap)**2: return False
        return True
    def add(self,x,r,M,kind):
        x=np.asarray(x,float); idx=len(self.X); self.X.append(x); self.r.append(float(r)); self.M.append(float(M)); self.kind.append(kind)
        self.grid.setdefault(self._key(x),[]).append(idx); self.rmax=max(self.rmax,float(r)); return idx
    def arrays(self):
        return (np.array(self.X).reshape(-1,3),np.array(self.r),np.array(self.M),np.array(self.kind))
    def copy(self):
        p=Packer(self.h); p.grid={k:list(v) for k,v in self.grid.items()}; p.X=list(self.X); p.r=list(self.r); p.M=list(self.M); p.kind=list(self.kind); p.rmax=self.rmax; return p

def check_no_overlap(X,r,gap=GAP):
    """brute-force verification (KD tree); returns min of d - Ri - Rj over all pairs (must be >= gap)."""
    from scipy.spatial import cKDTree
    if len(X)<2: return np.inf
    t=cKDTree(X); pairs=t.query_pairs(2*r.max()+gap,output_type="ndarray")
    if len(pairs)==0: return np.inf
    d=np.linalg.norm(X[pairs[:,0]]-X[pairs[:,1]],axis=1)-r[pairs[:,0]]-r[pairs[:,1]]
    return float(d.min())

def rand_dir(rng,n=None):
    u=rng.normal(size=(3,) if n is None else (n,3)); return u/np.linalg.norm(u,axis=-1,keepdims=True)

def rand_in_shell(rng,r0,r1,n):
    return rand_dir(rng,n)*(rng.uniform(r0**3,r1**3,n)**(1/3))[:,None]

def unit_pattern(X):
    """unit-mass (1 kg) pattern of each position; returns (n,360)."""
    X=np.atleast_2d(X); out=np.empty((len(X),360))
    for i in range(len(X)): out[i]=pattern(X[i:i+1],np.array([1.0]))
    return out

def unit_pattern_fast(X):
    """vectorised version through field3 (same formula, checked)."""
    import lattice_core as lc
    P=np.atleast_2d(X)/RV; Dd=np.linalg.norm(P,axis=1)
    out=np.empty((len(P),360))
    for i in range(0,len(P),3000):
        sl=slice(i,i+3000); out[sl]=(1.0/KG)*lc.field3(np.arctan2(P[sl,1],P[sl,0]),np.arcsin(np.clip(P[sl,2]/Dd[sl],-1,1)),Dd[sl]).reshape(-1,360)
    return out

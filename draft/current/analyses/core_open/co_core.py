"""CORE OPEN: the CoE interior (within 43.9 Gpc of the CoE) is fully open to bodies -- no wall, no smooth-density cap.
Kernels are the existing ones, unchanged (imported read-only):
  outside / no-mirror : ../solution_space/ss_core.Model  (FC surround compensation)      -> class CO
  all-open            : ../solution_space/all_open/ao_core.AO (AO=open, FLOC compensation) -> class COAO
  giant partners      : ../solution_space/extreme_partners/ep_core.EP (strict net zero)   -> class COEP
Only change: how far out a body may sit, and the CoE-interior exclusion is dropped.
  reach along a direction U for a body of radius Rc:  dmax(U) = max(RMAX_OLD, core_exit(U) - Rc)
  where core_exit(U) is where the ray from us leaves the CoE-interior sphere (|x-CoE| = 43.9 Gpc).
  Toward the CoE the core runs from 14.3 to 102.1 Gpc from us, so bodies can reach the whole core,
  its far side included. Away from the core the old 43 Gpc (60 Gpc for partners) reach is kept.
Mass bounds unchanged: M 1e49..2e52 kg, deep voids delta=-0.7, Rv 0.5..4 Gpc. The void deficit follows the local
exp-lock density (inside the core rho is 1.02-1.07 of mean, so no effective cap there either)."""
import os,sys,numpy as np
ME=os.path.dirname(os.path.abspath(__file__)); CS=os.path.dirname(ME); SS=os.path.join(CS,"solution_space")
os.environ["AO"]="open"
sys.path.insert(0,SS); sys.path.insert(0,os.path.join(SS,"all_open")); sys.path.insert(0,os.path.join(SS,"extreme_partners"))
_c=os.getcwd()
import ss_core
from ss_core import *
import ao_core, ep_core
os.chdir(_c)
RMAX_OLD=43.0
CN=float(np.linalg.norm(C_G)); CU=C_G/CN

def core_exit(U):
    b=U@C_G; disc=b*b-(CN**2-A_IN_G**2)
    return np.where(disc>0,b+np.sqrt(np.maximum(disc,0)),0.0)
def dmax_dir(U,Rc,old=RMAX_OLD): return np.maximum(old,core_exit(U)-Rc)
def core_depth(X):
    """>0: inside the CoE interior by this many Gpc"""
    return A_IN_G-np.linalg.norm(np.atleast_2d(X)-C_G,axis=1)
def in_core(X): return core_depth(X)>0

def _dirs(lon,lat): return np.c_[np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]

class CO(Model):
    """outside our sphere, anywhere out to dmax(U) including the full CoE interior; no wall."""
    def geo(self,p):
        p=p.reshape(self.n,4); u,lon,lat,s=p.T; isC=self.k=="C"
        M=np.zeros(self.n); Rv=np.zeros(self.n); Rc=np.zeros(self.n)
        M[isC]=10**(np.log10(MLO)+s[isC]*np.log10(MHI/MLO)); Rv[~isC]=RLO+s[~isC]*(RHI-RLO); Rc[~isC]=Rv[~isC]
        Rc[isC]=(3*M[isC]/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3)
        U=_dirs(lon,lat); dmin=RMIN+Rc; dmx=dmax_dir(U,Rc)
        return u,U,M,Rv,Rc,dmin,dmx
    def decode(self,p):
        u,U,M,Rv,Rc,dmin,dmx=self.geo(p); X=U*(dmin+u*(dmx-dmin))[:,None]
        isV=self.k=="V"; M[isV]=-0.70*rho_gpc(X[isV])*vol(Rv[isV])*GPC**3
        return X,M,Rv,Rc
    def violations(self,p):
        X,M,Rv,Rc=self.decode(p); v=[np.zeros(1)]
        if self.n>1:
            i,j=np.triu_indices(self.n,1); v.append(np.maximum(0,Rc[i]+Rc[j]+GAP-np.linalg.norm(X[i]-X[j],axis=1)))
        return np.concatenate(v)
    def encode(self,X,M,Rv):
        return _encode(self,X,M,Rv,lambda Rc:RMIN+Rc)

class COAO(ao_core.AO):
    """all open (inside our sphere, shell, core); no density limit; ring-straddle rule as in all_open."""
    def decode(self,p):
        p=p.reshape(self.n,4); u,lon,lat,s=p.T; isC=self.k=="C"
        M=np.zeros(self.n); Rv=np.zeros(self.n); Rc=np.zeros(self.n)
        M[isC]=10**(np.log10(MLO)+s[isC]*np.log10(MHI/MLO)); Rv[~isC]=RLO+s[~isC]*(RHI-RLO); Rc[~isC]=Rv[~isC]
        Rc[isC]=(3*M[isC]/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3)
        U=_dirs(lon,lat); dmin=Rc+ao_core.GAP_US; dmx=dmax_dir(U,Rc); X=U*(dmin+u*(dmx-dmin))[:,None]
        isV=self.k=="V"; M[isV]=-0.70*rho_gpc(X[isV])*vol(Rv[isV])*GPC**3
        return X,M,Rv,Rc
    def encode(self,X,M,Rv):
        return _encode(self,X,M,Rv,lambda Rc:Rc+ao_core.GAP_US)

def _encode(m,X,M,Rv,dminf):
    k=m.k; X=np.array(X,float); M=np.array(M,float); Rv=np.array(Rv,float); p=np.zeros((len(k),4)); isC=k=="C"
    d=np.linalg.norm(X,axis=1); p[:,1]=np.arctan2(X[:,1],X[:,0]); p[:,2]=np.arcsin(np.clip(X[:,2]/d,-1,1))
    p[isC,3]=(np.log10(np.clip(np.abs(M[isC]),MLO,MHI))-np.log10(MLO))/np.log10(MHI/MLO); p[~isC,3]=(np.clip(Rv[~isC],RLO,RHI)-RLO)/(RHI-RLO)
    Rc=np.where(isC,(3*np.clip(np.abs(M),MLO,MHI)/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3),np.clip(Rv,RLO,RHI))
    U=X/d[:,None]; dmin=dminf(Rc); dmx=dmax_dir(U,Rc)
    p[:,0]=np.clip((d-dmin)/(dmx-dmin),0,1); return p.ravel()

class COEP(ep_core.EP):
    """giant partners + ordinary bodies; partners and bodies may sit in the CoE interior. Shell outer-edge rule for
    partners kept (|x-CoE|+Rp <= 72.5); partner reach max(60, core exit - Rp)."""
    def __init__(self,sign,npart,kinds=""):
        super().__init__(sign,npart,kinds)
        if kinds: self.ord=CO(kinds)
    def partners(self,p):
        _,q,_=self.split(p); u,lon,lat,s,dp=q.T; Rp=ep_core.RP_LO+s*(ep_core.RP_HI-ep_core.RP_LO)
        U=_dirs(lon,lat); dmin=RMIN+Rp; dmx=dmax_dir(U,Rp,ep_core.DMAX); X=U*(dmin+u*(dmx-dmin))[:,None]
        M=dp*rho_gpc(X)*vol(Rp)*GPC**3
        return X,M,Rp,dp
    def violations(self,p):
        du,Xp,Mp,Rp,dp,Xo,Mo,Rvo,Rco=self.decode(p); v=[]
        dc=np.linalg.norm(Xp-C_G,axis=1); v.append(np.maximum(0,dc+Rp-ep_core.SHELL_OUT))
        X=np.r_[Xp,Xo]; Rc=np.r_[Rp,Rco]
        if len(X)>1:
            i,j=np.triu_indices(len(X),1); v.append(np.maximum(0,Rc[i]+Rc[j]+GAP-np.linalg.norm(X[i]-X[j],axis=1)))
        return np.concatenate(v)

def random_core_point(rng,R=0.0):
    """uniform random point inside the CoE interior, at least R inside its edge and outside our sphere by R"""
    while True:
        x=C_G+rng.normal(size=3)*0; v=rng.normal(size=3); v/=np.linalg.norm(v); x=C_G+v*(A_IN_G-R)*rng.uniform()**(1/3)
        if np.linalg.norm(x)>=RMIN+R: return x

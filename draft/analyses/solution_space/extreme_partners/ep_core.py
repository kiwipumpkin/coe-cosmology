"""Extreme partners: our visible sphere (14.3 Gpc) gets a uniform density contrast delta_us (on the exp-lock density), and
1-3 giant partner regions outside it (in the shell, CoE interior excluded) carry the opposite (or free) sign.
  * sphere pattern: delta_us * S1 (sphere_unit.npz; only the density-gradient part survives, mostly m=1)
  * partner p: uniform sphere radius Rp in [1,15] Gpc, contrast delta_p (void: -0.9..0, overdensity: 0..+2),
    mass = delta_p * rho(centre) * V. Pull = point mass at its centre (exact outside it; ring is outside every partner).
    centre distance from us <= 60 Gpc; whole partner outside our sphere (|x| >= 14.3+Rp), outside the CoE interior
    (|x-CoE| >= 43.9+Rp) and inside the shell outer edge (|x-CoE|+Rp <= 72.5). Partners don't overlap anything.
  * optional ordinary objects (C masses / V deep voids) exactly as ../ss_core.Model.
  * NET MASS: strict -- sphere excess + partners + ordinary objects must sum to zero (|net| <= 1e51 kg, ~0.1% of the sphere mass).
    No smooth background compensation is used here: the partners themselves must do the compensating."""
import os,sys,numpy as np
ME=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.dirname(ME))
from ss_core import *
from scipy.optimize import least_squares
_d=np.load(os.path.join(ME,"sphere_unit.npz")); S1=_d["S1"]; MS1=float(_d["Ms1"])
RP_LO,RP_HI=1.0,15.0; DMAX=60.0; SHELL_OUT=72.5; NET_TOL=1e51
SIGN={"overwide":((0.0,0.6),(-0.9,0.0)),"over":((0.0,0.3),(-0.9,0.0)),"under":((-0.3,0.0),(0.0,2.0)),"free":((-0.3,0.3),(-0.9,2.0))}
class EP:
    def __init__(self,sign,npart,kinds=""):
        self.sign=sign; self.np_=npart; self.kinds=kinds; self.ord=Model(kinds) if kinds else None
        (dl,dh),(pl,ph)=SIGN[sign]; self.pl,self.ph=pl,ph
        lo=[dl]+[0,-4*np.pi,-np.pi/2,0,pl]*npart; hi=[dh]+[1,4*np.pi,np.pi/2,1,ph]*npart
        if kinds: lo+=list(self.ord.lo); hi+=list(self.ord.hi)
        self.lo=np.array(lo,float); self.hi=np.array(hi,float)
    def split(self,p):
        return p[0],p[1:1+5*self.np_].reshape(self.np_,5),p[1+5*self.np_:]
    def partners(self,p):
        _,q,_=self.split(p); u,lon,lat,s,dp=q.T; Rp=RP_LO+s*(RP_HI-RP_LO)
        dirs=np.c_[np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]
        dmin=RMIN+Rp; d=dmin+u*(DMAX-dmin); X=dirs*d[:,None]
        M=dp*rho_gpc(X)*vol(Rp)*GPC**3
        return X,M,Rp,dp
    def decode(self,p):
        du,_,po=self.split(p); Xp,Mp,Rp,dp=self.partners(p)
        if self.ord: Xo,Mo,Rvo,Rco=self.ord.decode(po)
        else: Xo=np.zeros((0,3)); Mo=np.zeros(0); Rvo=np.zeros(0); Rco=np.zeros(0)
        return du,Xp,Mp,Rp,dp,Xo,Mo,Rvo,Rco
    def pattern(self,p):
        du,Xp,Mp,Rp,dp,Xo,Mo,Rvo,Rco=self.decode(p)
        w=du*S1+Mp@unit_field(Xp)
        if len(Mo): w=w+Mo@unit_field(Xo)
        return w
    def net(self,p):
        du,Xp,Mp,Rp,dp,Xo,Mo,Rvo,Rco=self.decode(p); return du*MS1+Mp.sum()+Mo.sum()
    def violations(self,p):
        du,Xp,Mp,Rp,dp,Xo,Mo,Rvo,Rco=self.decode(p); v=[]
        dc=np.linalg.norm(Xp-C_G,axis=1)
        v.append(np.maximum(0,A_IN_G+Rp-dc)); v.append(np.maximum(0,dc+Rp-SHELL_OUT))
        X=np.r_[Xp,Xo]; Rc=np.r_[Rp,Rco]
        if len(X)>1:
            i,j=np.triu_indices(len(X),1); v.append(np.maximum(0,Rc[i]+Rc[j]+GAP-np.linalg.norm(X[i]-X[j],axis=1)))
        if self.ord: v.append(np.maximum(0,A_IN_G+Rco-np.linalg.norm(Xo-C_G,axis=1)))
        return np.concatenate(v)
    def leftover(self,p): return rms(R-self.pattern(p))/FLOOR
SC=1/(FLOOR*np.sqrt(360))
def fit(m,p0):
    lam=30.0; p=np.clip(p0,m.lo,m.hi)
    for it in range(5):
        f=lambda q: np.concatenate([(R-m.pattern(q))*SC,lam*m.violations(q),[lam*m.net(q)/NET_TOL*1e-3]])
        r=least_squares(f,p,bounds=(m.lo,m.hi),method="trf",x_scale="jac",max_nfev=80*len(p),xtol=1e-9,ftol=1e-9)
        p=r.x
        if m.violations(p).max(initial=0)<=1e-3 and abs(m.net(p))<=NET_TOL: break
        lam*=10
    return p

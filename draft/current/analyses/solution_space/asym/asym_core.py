"""Asymmetric (no-mirror) fits of masses + deep voids.  Reuses ../ss_core.py unchanged (read-only).
RULE (mirror = reflection across the plane through us, the CoE and the galactic pole, i.e. l -> 2*15.16 - l;
because the residual ring cannot tell b from -b, the mirror image may also have its latitude flipped):
  * Two objects of the SAME kind are a forbidden mirror match if ALL of
        angle between object j and the mirror image of object i  < 20 deg
        distances from us differ by < 15 %
        sizes differ by < 25 %   (mass for masses, radius for voids)
  * No object may lie within 10 deg of the mirror plane (that object would be its own mirror image;
    for objects in the galactic plane this is 'within 10 deg of the l = 15/195 axis').
(Tolerances chosen so the old near-mirror fits are excluded: the old N=8 best had pairs only 16 deg off-mirror.
 Sensitivity variant ASYM_ANG=10 -> 10 deg / 5 deg.)
Implemented in the fit as a smooth penalty (8-norm of the three normalised mismatches must be >= 1.17, which implies
at least one mismatch >= 1 i.e. outside the box) and checked afterwards with the hard box rule above."""
import os, sys, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.dirname(HERE))
from ss_core import *
from scipy.optimize import least_squares
ANG_TOL=float(os.environ.get("ASYM_ANG",20)); D_TOL=0.15; S_TOL=0.25; PLANE_TOL=ANG_TOL/2; PNORM=8; TARGET=1.17
MODE=os.environ.get("ASYM_MODE","box")   # "pos": stricter variant -- mirror position alone (angle < tol) is forbidden, whatever the distance/size
A=np.radians(2*AXIS); CA,SA=np.cos(A),np.sin(A)
NRM=np.array([-np.sin(np.radians(AXIS)),np.cos(np.radians(AXIS)),0.0])   # normal of the mirror plane
def mirror_xy(U): return np.c_[CA*U[:,0]+SA*U[:,1],SA*U[:,0]-CA*U[:,1],U[:,2]]
def pair_metrics(k,X,M,Rv):
    """returns list of (i,j,ang_deg,dfrac,sfrac) for same-kind pairs and plane angles per object"""
    d=np.linalg.norm(X,axis=1); U=X/d[:,None]; Ua=U.copy(); Ua[:,2]=np.abs(Ua[:,2]); Mi=mirror_xy(Ua)
    size=np.where(k=="C",np.abs(M),Rv); out=[]
    for i in range(len(k)):
        for j in range(i+1,len(k)):
            if k[i]!=k[j]: continue
            ch=np.linalg.norm(Mi[i]-Ua[j]); ang=np.degrees(2*np.arcsin(np.clip(ch/2,0,1)))
            out.append((i,j,ang,abs(np.log(d[j]/d[i])),abs(np.log(size[j]/size[i]))))
    plane=np.degrees(np.arcsin(np.clip(np.abs(U@NRM),0,1)))
    return out,plane
def asym_viol(k,X,M,Rv):
    pm,plane=pair_metrics(k,X,M,Rv); v=[]
    for i,j,ang,dl,sl in pm:
        if MODE=="pos": v.append(max(0.0,1.05-ang/ANG_TOL)); continue
        s=((ang/ANG_TOL)**PNORM+(dl/np.log(1+D_TOL))**PNORM+(sl/np.log(1+S_TOL))**PNORM)**(1/PNORM)
        v.append(max(0.0,TARGET-s))
    v+=list(np.maximum(0,(PLANE_TOL*1.05-plane)/PLANE_TOL))
    return np.array(v)
def hard_check(k,X,M,Rv):
    """hard box rule: returns (n_mirror_matches, n_on_plane, closest pair summary)"""
    pm,plane=pair_metrics(k,X,M,Rv)
    match=[(i,j) for i,j,a,dl,sl in pm if a<ANG_TOL and (MODE=="pos" or (dl<np.log(1+D_TOL) and sl<np.log(1+S_TOL)))]
    # closest-to-mirror pair by normalised max mismatch
    best=min(pm,key=lambda t:max(t[2]/ANG_TOL,t[3]/np.log(1+D_TOL),t[4]/np.log(1+S_TOL))) if pm else None
    return len(match),int((plane<PLANE_TOL).sum()),best,plane
SC=1/(FLOOR*np.sqrt(360))
def fit_one(m,p0,asym=True):
    lam=30.0; p=p0
    def allv(q):
        X,M,Rv,Rc=m.decode(q); v=m.violations(q)
        return v,(asym_viol(m.k,X,M,Rv) if asym else np.zeros(0))
    for it in range(5):
        f=lambda q: np.concatenate([(R-m.pattern(q))*SC, lam*m.violations(q), lam*(allv(q)[1])])
        r=least_squares(f,p,bounds=(m.lo,m.hi),method="trf",x_scale="jac",max_nfev=60*len(p),xtol=1e-8,ftol=1e-8)
        p=r.x; v,va=allv(p)
        if v.max(initial=0)<=1e-3 and va.max(initial=0)<=1e-3: break
        lam*=10
    return p
def odd_fraction(m,p):
    """share of the layout's own ring pattern that is antisymmetric about the axis (0 = perfectly mirror-symmetric)"""
    w=m.pattern(p); th=np.degrees(TH); ref=(2*AXIS-th)%360
    wm=np.interp(ref,th%360,w,period=360); odd=(w-wm)/2
    return float(rms(odd)/rms(w))

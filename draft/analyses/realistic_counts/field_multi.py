"""Minimum-norm smooth density field (1x target) over inside / outside / both regions, with optional one-sided (no-mirror) rule.
Reuses ../inside_visible/smooth_field/sf_core.py (same kernel, residual, Fourier data). Penalty = integral of (S^-1 delta)^2 dV."""
import os, sys, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); SF=os.path.join(os.path.dirname(HERE),"inside_visible","smooth_field")
sys.path.insert(0,SF); _c=os.getcwd()
from sf_core import *
from inside_core import C_G, A_IN_G
from ss_core import AXIS
from scipy.ndimage import gaussian_filter
os.chdir(_c)
NRM=np.array([-np.sin(np.radians(AXIS)),np.cos(np.radians(AXIS)),0.0])
class G:
    def __init__(self,H,rlo,rhi,half=0,plane_deg=10.0,coe=True):
        n=int(np.ceil(2*rhi/H))+3; n+=n%2; ax=(np.arange(n)-(n-1)/2)*H
        X,Y,Z=np.meshgrid(ax,ax,ax,indexing="ij"); P0=np.stack([X,Y,Z],-1); r=np.linalg.norm(P0,axis=-1)
        mask=(r>=rlo)&(r<=rhi)
        if coe: mask&=np.linalg.norm(P0-C_G,axis=-1)>=A_IN_G
        if half:
            sang=np.degrees(np.arcsin(np.clip((P0@NRM)/np.maximum(r,1e-9),-1,1)))
            mask&=(half*sang>=plane_deg)
        self.H=H; self.mask=mask; self.n=n; self.P=P0[mask]; self.vc=H**3
        self.m=rho_gpc(self.P)*(H*GPC)**3
        A=np.zeros((2*MM,len(self.P)))
        for i in range(0,len(self.P),20000):
            sl=slice(i,i+20000); A[:,sl]=(CS@unit_field(self.P[sl]).T)*self.m[sl]/FLOOR
        self.A=A
    def smooth(self,v,s):
        f=np.zeros(self.mask.shape); f[self.mask]=v
        return gaussian_filter(f,s/self.H,mode="constant",truncate=4.0)[self.mask]
def solve(grids,s,target=1.0,kappa=1e3):
    rows=[];mrow=[]
    for g in grids:
        rows.append(np.array([g.smooth(a,s) for a in g.A])/np.sqrt(g.vc)); mrow.append(g.smooth(g.m/1e50,s)/np.sqrt(g.vc))
    D=np.r_[np.concatenate(rows,1),np.concatenate(mrow)[None]*kappa]; K=D@D.T; t=np.r_[yR,0.0]
    def field(lam):
        u=np.linalg.solve(K+lam*np.eye(len(K)),t); v=D.T@u; out=[];i=0
        for g in grids:
            n=len(g.P); out.append(g.smooth(v[i:i+n],s)/np.sqrt(g.vc)); i+=n
        return out
    def left(ds): return float(np.sqrt(((yR-sum(g.A@d for g,d in zip(grids,ds)))**2).sum()+TAIL**2))
    lo,hi=-40,20; best=None
    if left(field(10**lo))>target: return None
    for _ in range(70):
        mid=(lo+hi)/2
        if left(field(10**mid))>target: hi=mid
        else: lo=mid
    ds=field(10**lo); return ds,left(ds)
def stats(grids,ds):
    out=[]
    for g,d in zip(grids,ds):
        mm=g.m*d; out.append(dict(rms=float(np.sqrt(np.mean(d**2))),peak=float(np.abs(d).max()),dmax=float(d.max()),dmin=float(d.min()),
            M_plus=float(mm[mm>0].sum()),M_minus=float(-mm[mm<0].sum()),
            V_over_gpc3=float((d>0.01).sum()*g.vc),V_under_gpc3=float((d<-0.01).sum()*g.vc),
            mean_d_over=float(d[d>0.01].mean()) if (d>0.01).any() else 0.0,mean_d_under=float(d[d<-0.01].mean()) if (d<-0.01).any() else 0.0,
            Vregion=float(len(d)*g.vc)))
    return out

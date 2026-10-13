"""Continuous density-contrast field delta(x) inside the visible sphere, scored with the same Newtonian tidal kernel
(ss_core.unit_field == lattice_core.field3) and the same residual R. Linear inverse problem solved in data space.
Grid: cartesian, spacing H (Gpc); allowed cells R_US <= |x| <= R_MAX. Mass per cell per unit delta = rho_exp-lock * H^3.
Data: Fourier coefficients m=1..MM of the ring pattern (cos,sin)/sqrt2 in floor units -> ||y|| = rms/floor (m=0 removed by the
kernel's mean subtraction); the tail m>MM of R is added in quadrature (it is ~0)."""
import os, sys, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.dirname(HERE)); _c=os.getcwd()
from inside_core import R,FLOOR,RV,rms,GPC,RHO,rho_gpc,TH,unit_field,amps,M_VIS
os.chdir(_c)
MM=20
CS=np.r_[[np.cos(m*TH) for m in range(1,MM+1)],[np.sin(m*TH) for m in range(1,MM+1)]]*(2/360)/np.sqrt(2)   # (40,360)
def four(v): return CS@v                       # rms(v) ~= ||four(v)|| (exact for band-limited, zero-mean v)
yR=four(R)/FLOOR; TAIL=np.sqrt(max(0,(rms(R)/FLOOR)**2-yR@yR))

class Grid:
    def __init__(self,H=0.3,r_us=0.5,r_max=13.8):
        n=int(np.ceil(2*(RV+1.5*0)/H))+1; n+=n%2; self.n=n; self.H=H
        ax=(np.arange(n)-(n-1)/2)*H; self.ax=ax
        X,Y,Z=np.meshgrid(ax,ax,ax,indexing="ij"); r=np.sqrt(X*X+Y*Y+Z*Z)
        self.mask=(r>=r_us)&(r<=r_max); self.r=r
        P=np.c_[X[self.mask],Y[self.mask],Z[self.mask]]; self.P=P
        self.m=rho_gpc(P)*(H*GPC)**3                          # kg per unit delta per cell
        A=np.zeros((2*MM,len(P)))
        for i in range(0,len(P),20000):
            sl=slice(i,i+20000); A[:,sl]=(CS@unit_field(P[sl]).T)*self.m[sl]/FLOOR
        self.A=A                                             # data per unit delta in each cell (floor units)
        self.vol=self.mask.sum()*H**3
    def tofield(self,v):
        f=np.zeros(self.mask.shape); f[self.mask]=v; return f
    def smooth(self,v,s):
        """Gaussian smoothing (sigma = s Gpc) of a masked-cell vector, result restricted to the mask (FFT, zero-padded)."""
        if s<=0: return v
        from scipy.ndimage import gaussian_filter
        return gaussian_filter(self.tofield(v),s/self.H,mode="constant",truncate=4.0)[self.mask]

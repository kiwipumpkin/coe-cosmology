"""Einasto-style periodic backgrounds through our visible sphere, scored with the same ring kernel as
../sequential_open/wave_bg.py (pull at ring point minus pull on us, along the sightline, mean removed). Read-only reuse.
Period: Einasto et al. 1997 (Nature 385, 139; astro-ph/9701018): lambda0 = 2pi/k0 = 120 +/- 15 h^-1 Mpc (k0=0.052 h/Mpc),
grid step 120 +/- 20 h^-1 Mpc.  h=0.7 -> 171.4 Mpc = 0.1714 Gpc.
Models
  chess   : 3D chessboard, cubic cells of edge lambda/2 = 85.7 Mpc alternating +A / -A (square wave, rms = A)
  einasto : void at each cell centre, clusters on walls/edges/corners, cell edge = lambda = 171.4 Mpc;
            delta = -A (cos kx + cos ky + cos kz)/sqrt(1.5)   (rms = A; void centre -2.45A, corners +2.45A)
Grid axes: supergalactic (as Einasto Fig 1): SGX->(l137.37,b0), SGY->(l227.37,b83.68), SGZ->(l47.37,b6.32); ring plane = x-y (angle = l).
"""
import os,sys,numpy as np,numba as nb
sys.path.insert(0,"/workspace/clump_scan/sequential_open"); _c=os.getcwd()
from common import *          # C_G,A_IN_G,R,FLOOR,rms,RMIN, seq_core
os.chdir(_c)
KG=seq_core.KG; RV=seq_core.RV; GPC=seq_core.GPC; RHO=seq_core.RHO
H0h=0.7; LAM=0.120/H0h        # Gpc
AXIS=15.16
def gal(l,b):
    l,b=np.radians(l),np.radians(b); return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
Q_SG=np.c_[gal(137.37,0),gal(227.37,83.68),gal(47.37,6.32)]   # columns: SG axes in galactic coords
def rand_Q(rng):
    q,_=np.linalg.qr(rng.normal(size=(3,3))); return q
def delta(P,model,A,Q=Q_SG,phase=(0,0,0)):
    X=P@Q+np.asarray(phase)*LAM          # lattice-frame coords (Gpc)
    if model=="chess":
        c=LAM/2; s=(np.floor(X[:,0]/c)+np.floor(X[:,1]/c)+np.floor(X[:,2]/c))%2
        return A*(1-2*s)
    k=2*np.pi/LAM; return -A*(np.cos(k*X[:,0])+np.cos(k*X[:,1])+np.cos(k*X[:,2]))/np.sqrt(1.5)

@nb.njit(parallel=True,fastmath=True)
def _kern(x,y,z,m,e2,cth,sth):
    nt=len(cth); out=np.zeros(nt)
    for j in nb.prange(nt):
        c=cth[j]; s=sth[j]; acc=0.0
        for i in range(len(x)):
            dx=x[i]-c; dy=y[i]-s; dz=z[i]
            r2=dx*dx+dy*dy+dz*dz+e2
            acc+=m[i]*(dx*c+dy*s)/(r2*np.sqrt(r2))
        out[j]=acc
    return out
def kern(P,M,eps,nth=360):
    th=np.radians(np.arange(nth)*360/nth); cth,sth=np.cos(th),np.sin(th)
    X=P/RV; m=M/KG; e2=(eps/RV)**2
    a=_kern(X[:,0].copy(),X[:,1].copy(),X[:,2].copy(),m,e2,cth,sth)
    D3=(np.sum(X*X,1)+e2)**1.5; g=(m/D3)@X                # pull on us (vector)
    out=a-(g[0]*cth+g[1]*sth)
    return out-out.mean(), -(g[0]*cth+g[1]*sth)

def points(DT=1.0,h=LAM/6,sphere_only=False,local=1.0):
    """sample points in the torus within DT of the ring, plus a ball of radius `local` around us"""
    n=int(np.ceil((RV+DT)/h)); g=(np.arange(-n,n)+0.5)*h; nz=int(np.ceil(DT/h)); gz=(np.arange(-nz,nz)+0.5)*h
    XX,YY=np.meshgrid(g,g,indexing="ij"); rc=np.hypot(XX,YY); out=[]
    for zi in gz:
        sel=(rc-RV)**2+zi**2<DT**2
        if sphere_only: sel&=rc**2+zi**2<RV**2
        out.append(np.c_[XX[sel],YY[sel],np.full(sel.sum(),zi)])
    gl=(np.arange(-int(np.ceil(local/h)),int(np.ceil(local/h)))+0.5)*h
    L=np.stack(np.meshgrid(gl,gl,gl,indexing="ij"),-1).reshape(-1,3); L=L[np.linalg.norm(L,axis=1)<local]
    return np.concatenate(out+[L]),h

def pattern(model,A,Q=Q_SG,phase=(0,0,0),DT=1.0,h=None,nth=360,sphere_only=False,pts=None):
    if pts is None: P,h=points(DT,h or LAM/6,sphere_only)
    else: P,h=pts
    M=delta(P,model,A,Q,phase)*seq_core.rho_gpc(P)*GPC**3*h**3
    return kern(P,M,h,nth)

def fourier(w):
    n=len(w); th=np.radians(np.arange(n)*360/n)-np.radians(AXIS); m=np.arange(1,n//2)
    a=2*np.mean(w[None]*np.cos(m[:,None]*th),1); b=2*np.mean(w[None]*np.sin(m[:,None]*th),1)
    return m,a,b       # a: mirror-symmetric about the CoE axis, b: antisymmetric
def band(w,lo,hi):
    m,a,b=fourier(w); s=(m>=lo)&(m<=hi); return np.sqrt(np.sum(a[s]**2+b[s]**2)/2)/FLOOR
def sym_share(w,lo=1,hi=10**6):
    m,a,b=fourier(w); s=(m>=lo)&(m<=hi); return float(np.sum(a[s]**2)/np.sum(a[s]**2+b[s]**2))
def lowpass(w,mmax=12):
    m,a,b=fourier(w); n=len(w); th=np.radians(np.arange(n)*360/n)-np.radians(AXIS); s=m<=mmax
    return (a[s,None]*np.cos(m[s,None]*th)+b[s,None]*np.sin(m[s,None]*th)).sum(0)

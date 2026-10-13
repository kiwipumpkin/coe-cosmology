"""Smooth-field backgrounds on the filled region (core + our sphere + shell = 'both'):
  wave    : concentric shells about the CoE, delta = A cos(2 pi r_C / lam + phi), lam = 0.15 Gpc (BAO-like), with jitter
  lattice : jittered quasi-cubic lattice, delta = A/3 (cos kx' + cos ky' + cos kz'), spacing 0.12 Gpc (Einasto 1997), x' = jittered coords
Amplitude: peak contrast A = 0.15 (rms of matter contrast on ~60-75 Mpc scales, i.e. the realistic web smoothed to the wave/lattice scale);
the pattern is linear in A, so other amplitudes scale directly.
Jitter: smooth random phase / displacement field (Fourier modes with 0.5-2 Gpc wavelengths), rms 0.25 of a wavelength.
How the pull is computed:
  wave, coherent part : exact (shell theorem) for the spherical wave everywhere, delta M_enc(r) analytic
  near the ring       : the jittered field minus the unjittered field, summed cell by cell (0.0375 Gpc cells, Plummer softening = 1 cell)
                        within 1.0 Gpc of the ring (the high-m texture is set there; jitter farther out averages away)
  lattice             : periodic, so its far field cancels; summed cell by cell within 1.0 Gpc of the ring
Kernel: same tidal formula as field3 (pull at ring point minus pull on us, along the sightline, mean removed).
usage: python3 wave_bg.py wave|lattice SEED -> final_<kind>_s<seed>.npz (w, plus a ring-plane slab of the field for the side views)"""
import os,sys,time,numpy as np
from common import *
KG=seq_core.KG; RV=seq_core.RV; GPC=seq_core.GPC; RHO=seq_core.RHO
kind=sys.argv[1]; seed=int(sys.argv[2]); rng=np.random.default_rng(seed)
A=float(os.environ.get("AMP",0.15)); LAM=0.15 if kind=="wave" else 0.12; k=2*np.pi/LAM; H=0.0375; DT=1.0; JIT=0.25
TH=np.radians(np.arange(360)); cth,sth=np.cos(TH),np.sin(TH); RHO_G=RHO*GPC**3   # kg per Gpc^3
# smooth jitter field: sum of 40 random plane waves, wavelengths 0.5-2 Gpc
NW=40; kv=rng.normal(size=(NW,3)); kv/=np.linalg.norm(kv,axis=1,keepdims=True); kv*=(2*np.pi/rng.uniform(0.5,2.0,NW))[:,None]; ph=rng.uniform(0,2*np.pi,NW)
def jit(P,c=0):   # unit-rms smooth random field (3 independent fields for c=0,1,2)
    return np.sqrt(2/NW)*np.cos(P@np.roll(kv,c,0).T+np.roll(ph,7*c)).sum(1)
Q,_=np.linalg.qr(rng.normal(size=(3,3)))     # random lattice orientation
def delta(P,jitter=True):
    if kind=="wave":
        r=np.linalg.norm(P-C_G,axis=1); return A*np.cos(k*r+(2*np.pi*JIT*jit(P) if jitter else 0.0))
    Pp=P@Q+(LAM*JIT*np.c_[jit(P,0),jit(P,1),jit(P,2)] if jitter else 0.0)
    return A/3*(np.cos(k*Pp[:,0])+np.cos(k*Pp[:,1])+np.cos(k*Pp[:,2]))
def kern(P,M,eps):
    """(360,) pattern of point masses M (kg) at P (Gpc), softened; units as R (field3/KG convention)"""
    X=P/RV; e2=(eps/RV)**2; out=np.zeros(360)
    for i in range(0,len(X),20000):
        x=X[i:i+20000]; m=M[i:i+20000]/KG
        dx=x[:,0:1]-cth[None]; dy=x[:,1:2]-sth[None]; dz=x[:,2:3]
        r3=(dx*dx+dy*dy+dz*dz+e2)**1.5; D3=(np.sum(x*x,1)+e2)**1.5
        out+=m@((dx*cth[None]+dy*sth[None])/r3)-((m/D3)[:,None]*(x[:,0:1]*cth[None]+x[:,1:2]*sth[None])).sum(0)
    return out-out.mean()
t0=time.time()
# --- coherent spherical part (wave only), exact
w=np.zeros(360)
if kind=="wave":
    def Menc(r): return 4*np.pi*RHO_G*A*(r*r*np.sin(k*r)/k+2*r*np.cos(k*r)/k**2-2*np.sin(k*r)/k**3)
    ring=np.c_[cth,sth,np.zeros(360)]*RV; rt=np.linalg.norm(ring-C_G,axis=1); ru=np.linalg.norm(C_G)
    X=C_G/RV; dx=X[0]-cth; dy=X[1]-sth; dz=X[2]
    # pull toward C at ring point = Menc(r_t) * (C - p)/|C-p|^3 ; at us = Menc(|C|) * C/|C|^3  (same as field3 with M -> Menc at each end)
    ar=(Menc(rt)/KG)*(dx*cth+dy*sth)/((dx*dx+dy*dy+dz*dz)**1.5)-(Menc(ru)/KG)*(X[0]*cth+X[1]*sth)/np.sum(X*X)**1.5
    w+=ar-ar.mean()
w_coh=w.copy()
# --- near-ring cells: torus within DT of the ring, grid of H
n=int(np.ceil((RV+DT)/H)); g=(np.arange(-n,n+1)+0.5)*H; nz=int(np.ceil(DT/H)); gz=(np.arange(-nz,nz)+0.5)*H
w_near=np.zeros(360); ncell=0
for zi in gz:
    XX,YY=np.meshgrid(g,g,indexing="ij"); rho_c=np.hypot(XX,YY); sel=(rho_c-RV)**2+zi**2<DT**2
    P=np.c_[XX[sel],YY[sel],np.full(sel.sum(),zi)]
    d1=delta(P,True); d0=delta(P,False) if kind=="wave" else 0.0
    M=(d1-d0)*seq_core.rho_gpc(P)*GPC**3*H**3; w_near+=kern(P,M,H); ncell+=len(P)
w=w_coh+w_near
# slab of the field in the ring plane for the side views (zoom 2 Gpc box around the ring point at l=195, and a coarse full view)
zx=np.linspace(-1.5,1.5,301); ZX,ZY=np.meshgrid(zx,zx,indexing="ij"); c0=np.array([np.cos(np.radians(195)),np.sin(np.radians(195)),0])*RV
Pz=np.c_[ZX.ravel()+c0[0],ZY.ravel()+c0[1],np.zeros(ZX.size)]; zoom=delta(Pz,True).reshape(ZX.shape)
np.savez(f"{OUT}/final_{kind}_s{seed}.npz",w=w,w_coh=w_coh,w_near=w_near,zoom=zoom,zoom_center=c0,zoom_half=1.5,A=A,lam=LAM,ncell=ncell)
print(f"DONE {kind} s{seed} A={A} lam={LAM} cells={ncell:,} left={rms(R-w)/FLOOR:.2f} own={rms(w)/FLOOR:.2f} coherent={rms(w_coh)/FLOOR:.2f} near={rms(w_near)/FLOOR:.2f} t={time.time()-t0:.0f}s",flush=True)

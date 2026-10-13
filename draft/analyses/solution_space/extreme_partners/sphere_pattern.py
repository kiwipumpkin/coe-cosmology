"""Ring pattern of a uniform density contrast delta_us=1 filling our 14.3 Gpc visible sphere (density = exp-lock rho(x)).
Uniform-density part gives exactly zero tidal pattern on the ring (spherical shell theorem; pattern is mean-subtracted);
only the gradient part rho(x)-rho(0) contributes, integrated numerically. -> sphere_unit.npz (S1 pattern, mass Ms1)"""
import sys,os,numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","..","local_web"))
from web_core import pattern,rho_gpc,RMIN,GPC,FLOOR,rms,R
def grid(h):
    g=np.arange(-RMIN+h/2,RMIN,h); X,Y,Z=np.meshgrid(g,g,g,indexing="ij"); P=np.c_[X.ravel(),Y.ravel(),Z.ravel()]
    return P[np.linalg.norm(P,axis=1)<RMIN-1e-9]
out={}
for h in (0.5,0.35,0.25):
    P=grid(h); dV=h**3
    rho=rho_gpc(P); rho0=rho_gpc(np.zeros((1,3)))[0]
    w_grad=pattern(P,(rho-rho0)*dV*GPC**3); w_uni=pattern(P,np.full(len(P),rho0*dV*GPC**3))
    Ms=(rho*dV*GPC**3).sum()
    print(f"h={h} n={len(P)} M={Ms:.3e} kg  grad-part rms={rms(w_grad)/FLOOR:.4f}x floor  uniform-part (should be 0) rms={rms(w_uni)/FLOOR:.4f}x",flush=True)
    out[h]=(w_grad,Ms)
w,Ms=out[0.25]
print("rel change 0.35->0.25:",rms(out[0.35][0]-w)/rms(w))
print("exact volume mass rho0*V:",rho_gpc(np.zeros((1,3)))[0]*4/3*np.pi*RMIN**3*GPC**3)
np.savez("sphere_unit.npz",S1=w,Ms1=Ms)

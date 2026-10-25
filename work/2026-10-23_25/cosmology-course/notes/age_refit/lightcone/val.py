# validation: homogeneous FLRW limit (M->0, uniform density, p=1 -> E ∝ r0^2, simultaneous bang = open FLRW).
import numpy as np; from lc import *
from scipy.integrate import quad
for Om in (0.3,1e-4):
  # a(t): R=r0 a(t); E=0.5 V0^2 r0^2 ; m= (4pi/3) rho r0^3 G/c^2 -> H0^2 Om/2 r0^3
  H0=1.0; V0=np.sqrt(1-Om)*H0; mu=0.5*Om*H0**2
  P=dict(kind='pl',M=1e-30,rhoa=mu*3/(4*np.pi)*F.C**2/F.GU,L=1e9,s=1e-9,V0=V0,p=1.0)
  L=LC(P)
  # find tn where a=1: R(tn,0)=1
  from scipy.optimize import brentq
  tn=brentq(lambda t: L.R(t,np.array([0.]))[0]-1,1e-3,10)
  o=obs(L,tn); 
  for th in (0,90,180):
    r=ray(L,tn,np.radians(th),zmax=3,rtol=1e-10); z=np.expm1(r['lnz1']); lb=tn-r['t']
    # analytic lookback time t(z)
    tl=lambda zz: quad(lambda u: 1/((1+u)*np.sqrt(Om*(1+u)**3+(1-Om)*(1+u)**2)),0,zz)[0]
    err=max(abs(np.interp(zz,z,lb)-tl(zz))/tl(zz) for zz in (0.1,0.5,1,2.5))
    print(f'Om={Om} H={o["Hbar"]:.6f} (exp 1) th={th} max rel err lookback(z)={err:.2e} stop={r["stop"]}')

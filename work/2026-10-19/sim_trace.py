# Simultaneous recombination at 4c, no mass: same v11 light tracing, but stop each backward ray at emission time t_emit = our crossing time (instead of at R_rec).
import sys,numpy as np,base as V
u=4.0;V.uD=u
import v11
from scipy.integrate import solve_ivp
# optional args: core mass fraction of critical, output file. Default (no args) reproduces the original no-mass run.
frac=float(sys.argv[1]) if len(sys.argv)>1 else 0.0;OUT=sys.argv[2] if len(sys.argv)>2 else 'z4_sim.npz'
from scipy.optimize import brentq
M=frac*0.5*u**2;R=0.5;F=v11.field(M,R);wo,t0,T,W,S=F
te=R/u if M==0 else brentq(lambda t:V.Rfun(t,wo,M)-R,T[0]*1.01,t0)   # our own crossing time of R_rec (same as zsave.py)
def AB(t,w):
    Rr=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0]
    return Rr,Rp*Rp,2*Rp*Rtp,2*Rp*Rpp,Rr*Rr,2*Rr*Rt,2*Rr*Rp
def rhs(l,y):
    t,w,ph,pt,pw,pp=y;Rr,A,At,Aw,B,Bt,Bw=AB(t,w)
    return [-pt,pw/A,pp/B,0.5*At/A**2*pw**2+0.5*Bt/B**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B**2*pp**2,0]
hit=lambda l,y:y[0]-te;hit.terminal=True
A0=AB(t0,wo)[1];B0=AB(t0,wo)[4];out=[]
for ps in np.arange(0,181,1.0):
    pr=np.radians(ps);y0=[t0,wo,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sqrt(B0)*np.sin(pr)]
    s=solve_ivp(rhs,(0,-200),y0,events=[hit],rtol=1e-9,atol=1e-11)
    y=s.y_events[0][0] if s.t_events[0].size else np.full(6,np.nan)
    re=V.Rfun(y[0],y[1],M);out.append((ps,-y[3],y[0],t0-y[0],y[1],y[2],re))
o=np.array(out);np.savez(OUT,o=o,t0=t0,te=te,wo=wo,M=M)
print('M',M,'t0',t0,'t0 Gyr',t0*42*3.26156,'t_emit Gyr',te*42*3.26156,'t_emit',te,'resolved',np.isfinite(o[:,1]).sum(),'1+z min/max',np.nanmin(o[:,1]),np.nanmax(o[:,1]),'t0/te',t0/te)
G=42;print('emission radius from CoE (Gpc) min/max',np.nanmin(o[:,6])*G,np.nanmax(o[:,6])*G,'angle from CoE axis max deg',np.degrees(np.nanmax(np.abs(o[:,5]))))
xe=o[:,6]*np.cos(o[:,5]);ye=o[:,6]*np.sin(o[:,5]);print('distance emission point to our crossing point (Gpc) min/max',np.nanmin(np.hypot(xe-V.Rfun(te,wo,M),ye))*G,np.nanmax(np.hypot(xe-V.Rfun(te,wo,M),ye))*G)

# Copy of famP.py model build (unchanged physics), but instead of tracing to recombination it traces
# low-z null geodesics from our position and returns the angular-diameter distance D_A(psi,z) for z<=0.47.
# psi=0 is the ray looking AWAY from the CoE (backward ray moves to larger x); CoE is at psi=180.
# usage: famZ.py fam p frac out   (env as famP: SL, C1..C4, XS, DRV)
import sys,os,numpy as np
src=open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/constrained/famP.py').read()
exec(src.split('def rho(')[0])
from scipy.integrate import solve_ivp
def AB(t,w):
    Rr=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0]
    return Rr,Rp*Rp,2*Rp*Rtp,2*Rp*Rpp,Rr*Rr,2*Rr*Rt,2*Rr*Rp
def rhs(l,y):
    t,w,ph,pt,pw,pp=y;Rr,A,At,Aw,B_,Bt,Bw=AB(t,w)
    return [-pt,pw/A,pp/B_,0.5*At/A**2*pw**2+0.5*Bt/B_**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B_**2*pp**2,0]
ZMAX=0.47;ZG=np.geomspace(0.003,ZMAX,40);PS=np.linspace(1,179,45);dps=np.radians(0.02)
ev=lambda l,y:-y[3]-1-ZMAX-0.01;ev.terminal=True
A0=AB(t0,1.0)[1];B0=AB(t0,1.0)[4]
def trace(pr):
    y0=[t0,1.0,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sqrt(B0)*np.sin(pr)]
    s=solve_ivp(rhs,(0,-50),y0,events=[ev],rtol=1e-10,atol=1e-13,dense_output=True)
    ll=np.linspace(0,s.t[-1],3000);Y=s.sol(ll);zz=-Y[3]-1
    return np.array([np.interp(ZG,zz,Y[k]) for k in range(6)])
DA=np.zeros((len(PS),len(ZG)));XE=np.zeros_like(DA)
for i,ps in enumerate(PS):
    pr=np.radians(ps);a=trace(pr);b=trace(pr+dps)
    for j in range(len(ZG)):
        t,w,ph,pt,pw,pp=a[:,j];Rr,A,_,_,B_,_,_=AB(t,w)
        dw=b[1,j]-w;dph=b[2,j]-ph;kw=pw/A;kp=pp/B_;nk=np.sqrt(A*kw*kw+B_*kp*kp)
        xs2=A*dw*dw+B_*dph*dph;xn=(A*dw*kw+B_*dph*kp)/nk
        Dpsi=np.sqrt(max(xs2-xn*xn,0))/dps;Dphi=Rr*np.sin(ph)/np.sin(pr)
        DA[i,j]=np.sqrt(Dpsi*Dphi);XE[i,j]=w
np.savez(out,DA=DA,XE=XE,ZG=ZG,PS=PS,ratio=ratio,fam=fam,p=p,frac=frac,DR=DR,C=np.array([float(os.environ.get('C%d'%k,'0')) for k in (1,2,3,4)]))
print(out,'ratio',ratio,'DA(psi=1,179 @z=.075)',np.interp(.075,ZG,DA[0]),np.interp(.075,ZG,DA[-1]),flush=True)

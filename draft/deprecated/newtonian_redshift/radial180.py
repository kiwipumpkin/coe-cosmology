import numpy as np
exec(open('ourframe14.py').read().split('def run')[0])
def ray(D,rs,NT=200000,Rlist=np.arange(3,14.31,1.0).tolist()+[14.3]):
    GM=rs*GPC*c**2/2; umax=3*(D+18)*GPC/t0+3*np.sqrt(GM/(0.1*GPC)) if GM>0 else 3*(D+18)*GPC/t0
    u=np.geomspace(1e-7,1,1500)*umax
    uo=np.interp(D*GPC,shells(GM,u,t0)[0],u)
    t=t0; r=D*GPC; lnz=0.; dt=t0/NT; out={}; minr=r; rows=[]
    Rnow=shells(GM,u,t0)[0]
    def V(r,t):
        Rt,Vt=shells(GM,u,t); return np.interp(r,Rt,Vt),np.interp(r,Rt,np.gradient(Vt)/np.gradient(Rt)),Rt
    while t>dt and len(out)<len(Rlist):
        v1,g1,_=V(r,t); rm=r-0.5*dt*(v1+c); v2,g2,_=V(rm,t-.5*dt)
        r=r-dt*(v2+c); lnz+=dt*g2; t-=dt; minr=min(minr,r)
        if r<=0: rows.append(("CROSSED r=0 at t=%.4f Gyr"%(t/3.156e16),)); break
        _,_,Rt=V(r,t); lab=np.interp(r,Rt,u); dnow=(D*GPC-np.interp(lab,u,Rnow))/GPC
        for R in Rlist:
            if R not in out and dnow>=R:
                out[R]=(c/1e3*lnz/(R*1e3), np.expm1(lnz), t/3.156e16, r/GPC)
    return out,rows
import sys
for D in (20,40,60):
  for rs in (0.,1.,10.,100.):
    out,rows=ray(D,rs,NT=4000)
    print("D=%d rs=%g"%(D,rs),rows)
    for R in sorted(out):
        Hc,z,te,re=out[R]; print("   R=%5.1f  H180=%.1f  z=%.2f  t_emit=%.3f Gyr  r_emit(CoE)=%.4f Gpc"%(R,Hc,z,te,re))

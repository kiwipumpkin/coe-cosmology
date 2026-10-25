# Zero-mass LTB: R(r,t)=v(r)(t-tB(r)), label r = areal radius today (t1). v(r)=uD (r/D)^p. Today R(r,t1)=r => tB(r)=t1-r/v(r).
# Metric ds^2=-dt^2+R'^2/(1+v^2)dr^2+R^2 dphi^2 (exact GR, M=0, E=v^2/2). p=1 & common bang = Milne.
import numpy as np,json,sys
from scipy.integrate import solve_ivp
D=1.0;uD=1.05;T1=D/uD
class U:
    def __init__(s,p):s.p=p
    def v(s,r):return uD*(r/D)**s.p
    def vp(s,r):return s.p*s.v(r)/r
    def vpp(s,r):return s.p*(s.p-1)*s.v(r)/r**2
    def tB(s,r):return T1-r/s.v(r)
    def tBp(s,r):return -(1-s.p)/s.v(r)            # d/dr [r/v] = (1-p)/v
    def tBpp(s,r):return (1-s.p)*s.vp(s,r)/s.v(r)**2 if False else (1-s.p)*s.p/(s.v(r)*r)
    def f(s,r,t):
        v,vp,vpp=s.v(r),s.vp(r),s.vpp(r);tb,tbp,tbpp=s.tB(r),s.tBp(r),s.tBpp(r);d=t-tb
        R=v*d;Rp=vp*d-v*tbp;Rpr=vpp*d-2*vp*tbp-v*tbpp;A=Rp**2/(1+v*v)
        At=2*Rp*vp/(1+v*v);Ar=2*Rp*Rpr/(1+v*v)-Rp**2*2*v*vp/(1+v*v)**2
        return R,Rp,v,vp,A,At,Ar
def rhs(U_,J):
    def F(l,y):
        t,r,ph,pt,pr=y;R,Rp,Rt,Rpt,A,At,Ar=U_.f(r,t);pph=J/R**2
        return [pt,pr,pph,-0.5*At*pr**2-R*Rt*pph**2,(0.5*Ar*pr**2+R*Rp*pph**2-(At*pt+Ar*pr)*pr)/A]
    return F
def trace(U_,th,event,lmax=50):
    R,Rp,v,vp,A,At,Ar=U_.f(D,T1);n=np.array([np.cos(th),np.sin(th)])
    pr=n[0]/np.sqrt(A);pph=n[1]/R;J=R*R*pph
    y0=[T1,D,0.0,-1.0,pr]
    ev=lambda l,y:event(U_,y);ev.terminal=True
    bad=lambda l,y:(y[0]-U_.tB(y[1]))-1e-9;bad.terminal=True
    small=lambda l,y:y[1]-1e-3;small.terminal=True
    so=solve_ivp(rhs(U_,J),(0,lmax),y0,method='DOP853',rtol=1e-10,atol=1e-12,events=[ev,bad,small],dense_output=True)
    if len(so.t_events[0])==0:return None
    y=so.y_events[0][0];l=so.t_events[0][0];t,r,ph,pt,_=y;R,Rp,v,vp,*_=U_.f(r,t)
    X=R*np.cos(ph);Y=R*np.sin(ph);Ob=np.array([D,0.]);L=so.sol(np.linspace(0,l,400))
    Rs=np.array([U_.f(rr,tt)[0] for tt,rr in zip(L[0],L[1])]);P=np.c_[Rs*np.cos(L[2]),Rs*np.sin(L[2])]
    return dict(opz=pt/-1.0,te=t,re=r,Re=R,theta_loc=(2*v/R+vp/Rp),travel=T1-t,Xe=X,Ye=Y,
                dist_emit_now=float(np.hypot(X-D,Y)),today_dist_of_emitter=float(np.hypot(r*np.cos(ph)-D,r*np.sin(ph))),path_areal_mixed=float(np.sum(np.hypot(*np.diff(P,axis=0).T))))
def ev_temp(arec):
    return lambda U_,y:(lambda R,Rp:np.cbrt(max(R,1e-300)**2*Rp/y[1]**2))(*U_.f(y[1],y[0])[:2])-arec
def ev_front(tmid,kap):   # emission time t = tmid + kap*(D - X): older (earlier) on anti-CoE side (X>D)
    return lambda U_,y:y[0]-(tmid+kap*(D-U_.f(y[1],y[0])[0]*np.cos(y[2])))
def lowz_H(U_,th,dl=1e-4):
    R,Rp,v,vp,A,At,Ar=U_.f(D,T1);n=np.array([np.cos(th),np.sin(th)]);pr=n[0]/np.sqrt(A);J=R*R*n[1]/R
    so=solve_ivp(rhs(U_,J),(0,dl),[T1,D,0,-1,pr],method='DOP853',rtol=1e-12,atol=1e-14);pt=so.y[3,-1]
    return (-pt-1)/dl    # dz/d(proper distance) at observer = local H along direction

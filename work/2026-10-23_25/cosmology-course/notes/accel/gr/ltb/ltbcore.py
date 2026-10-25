# LTB (spherical dust, exact GR), equatorial plane: ds^2 = -dt^2 + R'^2/(1+2E) dr^2 + R^2 dphi^2
# Units c=G=1, ring radius eps=1. r labels shells by their areal radius at the start time t1 (R(r,t1)=r).
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq, fsolve
from scipy.interpolate import RectBivariateSpline, CubicSpline
from scipy.interpolate import bisplev
def _ev(S,r,t):
    tx,ty,c=S.tck; kx,ky=S.degrees; return bisplev(r,t,(tx,ty,c,kx,ky))
class LTB:
    def __init__(s,D=10.0,Rdot=3.0,core_frac=0.0946,rmin=3.3,rmax=60.0,tspan=400.0,nr=500,nt=700):
        H=Rdot/D; Mcrit=H*H*D**3/2; s.Mc=core_frac*Mcrit; s.rhob=(Mcrit-s.Mc)/(4*np.pi/3*D**3)
        s.D=D; M=lambda r: s.Mc+4*np.pi/3*s.rhob*r**3; s.Mf=M
        ED=(Rdot**2-2*M(D)/D)/2; s.t1=s.age(D,M(D),ED)
        rg=np.linspace(rmin,rmax,nr); Eg=np.array([s.solveE(r,M(r)) for r in rg]); s.Es=CubicSpline(rg,Eg)
        tg=s.t1+np.r_[np.linspace(-0.5,0,6)[:-1], np.geomspace(1e-3,tspan,nt-6)]; tg=np.r_[s.t1-0.5,tg[1:]] if False else tg
        tg=np.unique(np.r_[s.t1+np.linspace(-0.5,0,6),tg]); Rg=np.empty((nr,len(tg)))
        for i,(r,E) in enumerate(zip(rg,Eg)):
            Mr=M(r); f=lambda t,R:[np.sqrt(max(2*Mr/R[0]+2*E,1e-30))]
            so=solve_ivp(f,(s.t1,tg[-1]),[r],t_eval=tg[tg>=s.t1],rtol=1e-12,atol=1e-12)
            sb=solve_ivp(f,(s.t1,tg[0]),[r],t_eval=tg[tg<s.t1][::-1],rtol=1e-12,atol=1e-12)
            Rg[i]=np.r_[sb.y[0][::-1],so.y[0]]
        s.Rs=RectBivariateSpline(rg,tg,Rg,kx=5,ky=5); s.rg=rg; s.tg=tg; s.Eg=Eg
        s.Rr=s.Rs.partial_derivative(1,0); s.Rt_=s.Rs.partial_derivative(0,1); s.Rrt=s.Rs.partial_derivative(1,1); s.Rrr=s.Rs.partial_derivative(2,0)
        s.minRp=s.Rs(rg,tg,dx=1).min()   # shell crossing check: R' must stay > 0
    def age(s,r,M,E): return quad(lambda x: r/np.sqrt(2*M/(r*x)+2*E),0,1,epsabs=1e-13,epsrel=1e-12,limit=400)[0]
    def solveE(s,r,M):
        lo=-M/r*(1-1e-9); return brentq(lambda E: s.age(r,M,E)-s.t1, lo, 1e6, xtol=1e-14)
    def fields(s,t,r):
        R=s.Rs.ev(r,t); Rp=_ev(s.Rr,r,t); Rt=_ev(s.Rt_,r,t); Rpt=_ev(s.Rrt,r,t); Rpp=_ev(s.Rrr,r,t)
        E=s.Es(r); Ep=s.Es(r,1); k=1+2*E
        G=Rp*Rp/k; Gt=2*Rp*Rpt/k; Gr=2*Rp*Rpp/k-Rp*Rp*2*Ep/k**2
        return R,Rp,Rt,G,Gt,Gr
    def rhs(s,t,y):
        r,ph,vr,vp,_=y; R,Rp,Rt,G,Gt,Gr=s.fields(t,r)
        Gtt=0.5*Gt*vr*vr+R*Rt*vp*vp
        ar=-(Gt/G*vr+0.5*Gr/G*vr*vr-R*Rp/G*vp*vp)+vr*Gtt
        ap=-(2*Rt/R*vp+2*Rp/R*vr*vp)+vp*Gtt
        return [vr,vp,ar,ap,-Gtt]
    def frame(s,t,r,ph,w):  # w=(wx,wy) local velocity in Cartesian-axis directions measured by comoving observer -> (vr,vphi)
        R,Rp,Rt,G,_,_=s.fields(t,r); wr=w[0]*np.cos(ph)+w[1]*np.sin(ph); wp=-w[0]*np.sin(ph)+w[1]*np.cos(ph)
        return np.array([wr/np.sqrt(G), wp/R])
    def kdu(s,t,r,vph,lnk,vp):
        R,Rp,Rt,G,_,_=s.fields(t,r); ut=1/np.sqrt(1-G*vp[0]**2-R*R*vp[1]**2)
        return np.exp(lnk)*ut*(-1+G*vph[0]*vp[0]+R*R*vph[1]*vp[1])
def cart2pol(X): return np.hypot(*X), np.arctan2(X[1],X[0])
def edge(L):
    ev=lambda t,y: min(y[0]-L.rg[0]-0.2, L.rg[-1]-1-y[0]); ev.terminal=True; return ev
def worldline(L,X0,w,ts,T):
    r,ph=cart2pol(X0); v=L.frame(ts,r,ph,w)
    so=solve_ivp(L.rhs,(ts,T),[r,ph,*v,0.0],method="DOP853",dense_output=True,rtol=1e-11,atol=1e-13,events=edge(L))
    tend=so.t[-1]
    return (lambda t: None if (t<ts-1e-12 or t>tend+1e-12) else so.sol(t)), so.status
def photon(L,ro,po,beta,T,tlow):
    R,Rp,Rt,G,_,_=L.fields(T,ro); v=[np.cos(beta)/np.sqrt(G), np.sin(beta)/R]   # beta from outward radial, comoving frame
    return solve_ivp(L.rhs,(T,tlow),[ro,po,*v,0.0],method="DOP853",dense_output=True,rtol=1e-11,atol=1e-13,events=edge(L))
def observe(L,ys,ems,T):
    ro,po,vro,vpo,_=ys; out=[]; prev=None
    Xc=lambda y: y[0]*np.array([np.cos(y[1]),np.sin(y[1])])     # comoving-label Cartesian (for the starting guess only)
    Xo=Xc(ys); aD=lambda t: L.Rs.ev(ro,t)/ro
    from scipy.integrate import quad as _q
    for em,ts in ems:
        if em(ts) is None: out.append(np.nan); continue
        hi=T-1e-6
        while em(hi) is None: hi=ts+(hi-ts)*0.99
        g=lambda s_: np.linalg.norm(Xo-Xc(em(s_)))-_q(lambda t: 1/aD(t),s_,T)[0]
        s0=brentq(g,ts,hi) if g(ts)<0<g(hi) else ts
        ye=em(s0); dX=Xc(ye)-Xo; b0=np.arctan2(-(dX@[-np.sin(po),np.cos(po)]),-(dX@[np.cos(po),np.sin(po)]))
        def res(z):
            b,s_=z
            if not (ts<=s_<=T): return [10+abs(s_),10]
            so=photon(L,ro,po,b,T,ts-0.2)
            if so.t[-1]>s_: return [5,5]
            e=em(s_)
            if e is None: return [7,7]
            y=so.sol(s_); Rr=L.Rs.ev(e[0],s_); return [y[0]-e[0], Rr*(y[1]-e[1])]
        sol=None
        for gs in ([b0,s0],prev):
            if gs is None: continue
            z=fsolve(res,gs,xtol=1e-12)
            if np.max(np.abs(res(z)))<1e-8: sol=z;break
        if sol is None: out.append(np.nan); continue
        prev=sol; b,s_=sol; so=photon(L,ro,po,b,T,ts-0.2); y=so.sol(s_); e=em(s_)
        R,Rp,Rt,G,_,_=L.fields(T,ro); vph=[np.cos(b)/np.sqrt(G), np.sin(b)/R]
        out.append(L.kdu(T,ro,vph,0.0,[vro,vpo])/L.kdu(s_,y[0],y[2:4],y[4],e[2:4]))
    return np.array(out)

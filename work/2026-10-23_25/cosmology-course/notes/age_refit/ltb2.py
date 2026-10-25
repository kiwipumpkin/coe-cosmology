"""Analytic-LTB + adaptive null geodesics (v2).  c=1 units (lengths Gpc, time Gpc/c).
Shells: R(t,x) from the exact parametric LTB/Kepler solution with smooth m(x) (CubicSpline of G M_enc/c^2 in x=ln r0),
E(x) and bang time tB(x) from a reference state (forward: t=0, R=r0, U=V0 (r0/r_i)^p; backward: t=t_n, R=r0, U=Uo (r0/D)^p).
Metric: ds^2=-dt^2+A^2 dr0^2+R^2 dphi^2, A=R'/W, W=sqrt(1+2E).  Derivatives R', Rdot', A' by centred differences of the
ANALYTIC solution in x (dx=1e-4) -> smooth (no table interpolation).
Rays: y=(x, phi, ln kt, q), q=A k^r/k^t (local radial direction cosine of propagation), k^phi=L/R^2, parameter t,
solve_ivp LSODA rtol=1e-10.  1+z = kt_emit/kt_obs (comoving dust emitters/observer)."""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp, cumulative_trapezoid
import forward as F
C=F.C
class LTB:
    def __init__(s,P,mode,r_ref_obs=1.0,t_ref=None):
        xs=np.linspace(np.log(1e-3),np.log(3e4),20001); r=np.exp(xs)
        rho=P['rhoa']*F.shape(P['kind'],r/P['L'],P.get('s'))
        mm=P['M']+np.r_[0,cumulative_trapezoid(4*np.pi*rho*r**3,xs)]      # kg (dm = 4 pi rho r^2 dr = 4 pi rho r^3 dx)
        s.m=CubicSpline(xs,F.GU*mm/C**2)                                      # Gpc
        s.mode=mode; s.P=P; s.xmin,s.xmax=xs[0],xs[-1]
        if mode=='fwd': s.Uref=lambda x: P['V0']/C*np.exp(x)**P['p']; s.tref=0.0
        else: s.Uref=lambda x: P['Uo']/C*(np.exp(x)/r_ref_obs)**P['p']; s.tref=t_ref
    def E(s,x): return 0.5*s.Uref(x)**2-s.m(x)/np.exp(x)
    def _T(s,x,R):   # time since bang to reach R on expanding branch
        mu=s.m(x); E=s.E(x)
        out=np.empty_like(R)
        rel=E*np.exp(x)/mu; pos=rel>1e-9; neg=rel<-1e-9; par=~(pos|neg)
        a=mu[pos]/(2*E[pos]); eta=np.arccosh(1+R[pos]/a); out[pos]=np.sqrt(a**3/mu[pos])*(np.sinh(eta)-eta)
        a=mu[neg]/(-2*E[neg]); c=np.clip(1-R[neg]/a,-1,1); eta=np.arccos(c)
        if s.mode=='back' or True:
            ur=s.Uref(x[neg]) if np.ndim(x) else s.Uref(x); eta=np.where(np.asarray(ur)<0,2*np.pi-eta,eta)
        out[neg]=np.sqrt(a**3/mu[neg])*(eta-np.sin(eta))
        out[par]=np.sqrt(2/(9*mu[par]))*R[par]**1.5
        return out
    def tB(s,x): return s.tref-s._T(x,np.exp(x))      # reference R=r0 at t_ref (expanding branch assumed)
    def R(s,t,x):
        x=np.atleast_1d(x).astype(float); t=np.broadcast_to(t,x.shape).astype(float)
        mu=s.m(x); E=s.E(x); tau=t-s.tB(x); out=np.full_like(x,np.nan); ok=tau>0
        rel=E*np.exp(x)/mu; pos=ok&(rel>1e-9); neg=ok&(rel<-1e-9); par=ok&~(pos|neg)
        if pos.any():
            a=mu[pos]/(2*E[pos]); T=tau[pos]/np.sqrt(a**3/mu[pos])      # sinh eta - eta = T
            eta=np.where(T<1,(6*T)**(1/3),np.arcsinh(T+1)+0.5)
            for _ in range(60): eta=eta-(np.sinh(eta)-eta-T)/(np.cosh(eta)-1+1e-300)
            out[pos]=a*(np.cosh(eta)-1)
        if neg.any():
            a=mu[neg]/(-2*E[neg]); T=tau[neg]/np.sqrt(a**3/mu[neg])      # full cycloid eta-sin(eta)=T, 0<T<2pi
            eta=np.where(T<1,(6*T)**(1/3),np.where(T>2*np.pi-1,2*np.pi-(6*np.clip(2*np.pi-T,0,None))**(1/3),T*0+np.pi))
            for _ in range(80): eta=np.clip(eta-(eta-np.sin(eta)-T)/(1-np.cos(eta)+1e-300),1e-12,2*np.pi-1e-12)
            Rn=a*(1-np.cos(eta)); Rn[T>=2*np.pi]=np.nan          # crunched: shell has fallen into the central mass
            out[neg]=Rn; s._eta_last=(neg,eta)
        if par.any(): out[par]=(4.5*mu[par])**(1/3)*tau[par]**(2/3)
        return out
    def Rdot(s,t,x):
        x=np.atleast_1d(x).astype(float); R=s.R(t,x); v=np.sqrt(np.clip(2*s.m(x)/R+2*s.E(x),0,None))
        E=s.E(x); neg=(E*np.exp(x)/s.m(x))<-1e-9
        if neg.any():   # collapsing half of the cycloid: tau beyond the turnaround time pi*sqrt(a^3/mu)
            a=s.m(x)/(-2*E+1e-300); tau=np.broadcast_to(t,x.shape)-s.tB(x); col=neg&(tau>np.pi*np.sqrt(np.abs(a)**3/s.m(x)))
            v=np.where(col,-v,v)
        return v
    def interior(s,x):  # E<-1/2: W^2=1+2E<0, no real exterior LTB slice -> inside the trapped (Schwarzschild-interior) region
        return 1+2*s.E(np.atleast_1d(x))<=0
    def rho_rel(s,t,x,dx=1e-4):
        x=np.atleast_1d(x).astype(float); Rp=(s.R(t,x+dx)-s.R(t,x-dx))/(2*dx)/np.exp(x); R=s.R(t,x)
        return s.m(x,1)/np.exp(x)/(R**2*Rp)        # proportional to rest-mass density (4 pi rho G/c^2 = m_r/(R^2 R'))
    def geo(s,t,x,dx=1e-4):
        R=s.R(t,x); Rp_=(s.R(t,x+dx)-s.R(t,x-dx))/(2*dx)/np.exp(x)
        U=s.Rdot(t,x); Up=(s.Rdot(t,x+dx)-s.Rdot(t,x-dx))/(2*dx)/np.exp(x)
        W=lambda xx: np.sqrt(1+2*s.E(xx))
        Af=lambda xx: (s.R(t,xx+dx)-s.R(t,xx-dx))/(2*dx)/np.exp(xx)/W(xx)
        A=Rp_/W(x); Ad=Up/W(x); Ap=(Af(x+dx)-Af(x-dx))/(2*dx)/np.exp(x)
        return R,U,Rp_,A,Ad,Ap
def ray(L,t0,x0,theta,t1,rtol=1e-10,atol=1e-12,dense=False):
    """Null ray observed at (t0,x0) from sky angle theta (0=away from CoE); integrate in t from t0 to t1 (t1<t0 backward).
    For forward use pass an initial state via ray_state()."""
    R,U,Rp,A,Ad,Ap=L.geo(t0,np.array([x0])); Lang=-R[0]*np.sin(theta)
    y0=[x0,0.0,0.0,-np.cos(theta)]
    return integ(L,Lang,t0,y0,t1,rtol,atol,dense)
def rhs(L,Lang):
    def f(t,y):
        x,phi,lkt,q=y
        if not(L.xmin+1e-3<x<L.xmax-1e-3): return [np.nan]*4
        R,U,Rp,A,Ad,Ap=[v[0] for v in L.geo(t,np.array([x]))]
        if not np.isfinite(R) or not np.isfinite(A) or A<=0: return [np.nan]*4
        kt=np.exp(lkt); kr=q*kt/A; kph=Lang/R**2
        dkt=(-A*Ad*kr**2-R*U*kph**2)/kt                    # dkt/dt
        dkr=(-2*(Ad/A)*kt*kr-(Ap/A)*kr**2+(R*Rp/A**2)*kph**2)/kt
        dx=kr/(kt*np.exp(x)); dA=Ad+Ap*np.exp(x)*dx
        dq=(dA*kr+A*dkr)/kt-q*dkt/kt
        return [dx,kph/kt,dkt/kt,dq]
    return f
def integ(L,Lang,t0,y0,t1,rtol=1e-10,atol=1e-12,dense=False):
    f=rhs(L,Lang)
    ev=lambda t,y: y[0]-(L.xmin+1.0); ev.terminal=True; ev.direction=-1
    y=np.array(y0,float); t=t0; ts=[];ys=[]; ncross=0
    for _ in range(4):
        s=solve_ivp(f,(t,t1),y,method='LSODA',rtol=rtol,atol=atol,events=ev)
        ts.append(s.t); ys.append(s.y)
        if s.status!=1: break
        # radial passage through the centre (r0 < e^{xmin+1}): continue on the opposite side, direction cosine flips
        t=s.t[-1]; y=s.y[:,-1].copy(); y[1]+=np.pi; y[3]=-y[3]; ncross+=1; t=t+(t1-t)*1e-12
    s.t=np.concatenate(ts); s.y=np.concatenate(ys,axis=1); s.ncross=ncross; return s

# Generic diagonal metric ds^2 = -A dt^2 + B (dx^2+dy^2), A,B functions of (t,rho), rho=|X| from mass at origin.
# McVittie (flat matter FLRW, a=(t/t0)^(2/3)) ; Schwarzschild-isotropic when H=0.
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, brentq
class Metric:
    def __init__(s,m,t0=None):  # t0=None -> static (a=1)
        s.m=m; s.t0=t0
    def aH(s,t):
        if s.t0 is None: return 1.0,0.0
        return (t/s.t0)**(2/3), 2/(3*t)
    def AB(s,t,X):
        a,H=s.aH(t); rho=np.hypot(*X); mu=s.m/(2*a*rho) if s.m else 0.0
        A=((1-mu)/(1+mu))**2; B=(1+mu)**4*a*a
        dAdmu=-4*(1-mu)/(1+mu)**3; dBdmu=4*(1+mu)**3*a*a
        rh=X/rho
        gA=dAdmu*(-mu/rho)*rh; At=dAdmu*(-mu*H)
        gB=dBdmu*(-mu/rho)*rh; Bt=dBdmu*(-mu*H)+2*H*B
        return A,B,gA,At,gB,Bt,mu
    def rhs(s,t,y):
        X=y[:2]; v=y[2:4]; A,B,gA,At,gB,Bt,mu=s.AB(t,X); vv=v@v
        Gi=gA/(2*B) + Bt*v/B + (2*v*(gB@v)-gB*vv)/(2*B)
        Gt=At/(2*A) + (gA@v)/A + Bt*vv/(2*A)
        acc=-Gi+v*Gt
        return [v[0],v[1],acc[0],acc[1],-Gt]
    def ut(s,t,X,v): A,B=s.AB(t,X)[:2]; return 1/np.sqrt(A-B*(v@v))
    def coordvel(s,t,X,w):  # w = velocity measured by local static/comoving observer
        A,B=s.AB(t,X)[:2]; return w*np.sqrt(A/B)
def worldline(Mt,X0,w0,t0,t1,rtol=1e-12):
    v0=Mt.coordvel(t0,X0,w0)
    ev=lambda t,y: Mt.AB(t,y[:2])[6]-0.9 if Mt.m else 1.0; ev.terminal=True
    so=solve_ivp(Mt.rhs,(t0,t1),[*X0,*v0,0.0],method="DOP853",dense_output=True,rtol=rtol,atol=rtol*1e-2,events=ev)
    tend=so.t[-1]
    def at(t):
        if t<t0-1e-12 or t>tend: return None
        y=so.sol(t); return y[:2],y[2:4]
    return at, so.status==1
def photon(Mt,Xo,beta,to,tlow,rtol=1e-12):
    A,B=Mt.AB(to,Xo)[:2]; v=np.sqrt(A/B)*np.array([np.cos(beta),np.sin(beta)])
    ev=lambda t,y: Mt.AB(t,y[:2])[6]-0.9 if Mt.m else 1.0; ev.terminal=True
    return solve_ivp(Mt.rhs,(to,tlow),[*Xo,*v,0.0],method="DOP853",dense_output=True,rtol=rtol,atol=rtol*1e-2,events=ev)
def kdotu(Mt,t,X,vph,lnkt,vp):
    A,B=Mt.AB(t,X)[:2]; return np.exp(lnkt)*Mt.ut(t,X,vp)*(-A+B*(vph@vp))
def observe(Mt,Xo,vo,ems,to,guesses=None):
    out=[];prev=None;info=[]
    for k,(em,ts) in enumerate(ems):
        g=lambda s: np.linalg.norm(Xo-em(s)[0])-(to-s)
        if em(ts) is None or g(ts)>0: out.append(np.nan); info.append(None); continue
        s0=brentq(g,ts,to-1e-6); d=em(s0)[0]-Xo; b0=np.arctan2(-d[1],-d[0])
        def res(z):
            b,s=z
            if not (ts<=s<=to): return [10+abs(s),10]
            so=photon(Mt,Xo,b,to,ts-0.5)
            if so.t[-1]>s: return [5,5]
            e=em(s)
            if e is None: return [7,7]
            return so.sol(s)[:2]-e[0]
        sol=None
        for gs in ([b0,s0],prev, None if guesses is None else guesses[k]):
            if gs is None: continue
            z=fsolve(res,gs,xtol=1e-12)
            if np.max(np.abs(res(z)))<1e-8: sol=z;break
        if sol is None: out.append(np.nan); info.append(None); continue
        prev=sol; b,s=sol; so=photon(Mt,Xo,b,to,ts-0.5); y=so.sol(s); Xe,ve=em(s)
        A,B=Mt.AB(to,Xo)[:2]; vph_o=np.sqrt(A/B)*np.array([np.cos(b),np.sin(b)])
        out.append(kdotu(Mt,to,Xo,vph_o,0.0,vo)/kdotu(Mt,s,y[:2],y[2:4],y[4],ve)); info.append((b,s))
    return np.array(out),info

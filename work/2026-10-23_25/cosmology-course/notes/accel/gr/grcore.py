import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve
th=np.radians(np.arange(0,361,1.0)); Mp=np.array([-2.0,0.0])
def polar(X): d=X-Mp; return np.hypot(*d), np.arctan2(d[1],d[0])
def cart(r,p): return Mp+r*np.array([np.cos(p),np.sin(p)])
def relsum(v0,u,e):  # bulk v0 along x plus expansion u along e (both in local static frame)
    G0=1/np.sqrt(1-v0*v0); up=u*e; den=1+v0*up[0]
    return np.array([(up[0]+v0)/den, up[1]/(G0*den)])
def fourvel(M,r,p,w):  # w: velocity (x,y components) in local static frame -> (u^t,u^r,u^phi)
    f=1-2*M/r; g=1/np.sqrt(1-w@w); vr=w[0]*np.cos(p)+w[1]*np.sin(p); vp=-w[0]*np.sin(p)+w[1]*np.cos(p)
    return np.array([g/np.sqrt(f), g*vr*np.sqrt(f), g*vp/r])
def static_vel(M,r,p,U):  # inverse: local static-frame velocity (x,y) from 4-velocity
    f=1-2*M/r; vr=U[1]/(np.sqrt(f)*U[0]*np.sqrt(f)); vp=U[2]*r/(U[0]*np.sqrt(f))
    return np.array([vr*np.cos(p)-vp*np.sin(p), vr*np.sin(p)+vp*np.cos(p)])
HOR=lambda M:(lambda t,y: y[0]-2*M*(1+1e-6)-1e-12)
def geodesic(M,X0,w0,t0,t1):
    """timelike geodesic from coordinate position X0, local-static velocity w0, at t0 -> t1. returns f(t)->(X,U)"""
    r0,p0=polar(X0); U=fourvel(M,r0,p0,w0); E=(1-2*M/r0)*U[0]; L=r0**2*U[2]
    def rhs(t,y):
        r=y[0]; f=1-2*M/r; s=f/E
        return [y[2]*s, L*s/r**2, s*(-M/r**2+L**2/r**3-3*M*L**2/r**4)]
    ev=HOR(M); ev.terminal=True
    so=solve_ivp(rhs,(t0,t1),[r0,p0,U[1]],dense_output=True,rtol=1e-11,atol=1e-13,events=ev if M>0 else None)
    tend=so.t[-1]
    def at(t):
        if (t-tend)*(t1-t0)>0: return None
        r,p,ur=so.sol(t); return cart(r,p), np.array([E/(1-2*M/r),ur,L/r**2])
    return at, (so.status==1)
def thrust(M,X0,w,t0,t1):
    """worldline held at constant velocity w relative to local static observers"""
    def rhs(t,y):
        r,p=y; f=1-2*M/r; vr=w[0]*np.cos(p)+w[1]*np.sin(p); vp=-w[0]*np.sin(p)+w[1]*np.cos(p)
        return [f*vr, np.sqrt(f)*vp/r]
    ev=HOR(M); ev.terminal=True
    so=solve_ivp(rhs,(t0,t1),list(polar(X0)),dense_output=True,rtol=1e-11,atol=1e-13,events=ev if M>0 else None)
    tend=so.t[-1]
    def at(t):
        if (t-tend)*(t1-t0)>0: return None
        r,p=so.sol(t); return cart(r,p), fourvel(M,r,p,w)
    return at, (so.status==1)
def ray(M,Xo,beta,to,tlow):
    """null geodesic arriving at Xo at time to; beta = propagation angle from r-hat in observer's static frame. Integrated backward."""
    r0,p0=polar(Xo); f0=1-2*M/r0; L=r0*np.sin(beta)/np.sqrt(f0)
    def rhs(t,y):
        r=y[0]; f=1-2*M/r
        return [f*y[2], f*L/r**2, f*(L**2/r**3-3*M*L**2/r**4)]
    ev=HOR(M); ev.terminal=True
    so=solve_ivp(rhs,(to,tlow),[r0,p0,np.cos(beta)],dense_output=True,rtol=1e-11,atol=1e-13,events=ev if M>0 else None)
    return so,L
def kdotu(M,r,pr,L,U): return -U[0]+pr/(1-2*M/r)*U[1]+L*U[2]
def observe(M,Xo,Uo,emitters,to,tlow):
    """emitters: list of worldline functions t->(X,U). returns T_o/T_e per emitter, and info"""
    ro,po=polar(Xo); out=[]; prev=None; bad=0
    for k,em in enumerate(emitters):
        g=lambda s: np.linalg.norm(Xo-em(s)[0])-(to-s) if em(s) is not None else 1.0
        try: s0=brentq(g,tlow+1e-9,to-1e-6)
        except Exception: s0=prev[1] if prev is not None else to-1
        Xe=em(s0)[0]; d=Xe-Xo; dr=d@np.array([np.cos(po),np.sin(po)]); dp=d@np.array([-np.sin(po),np.cos(po)])
        fo=1-2*M/ro; b0=np.arctan2(-dp,-dr/np.sqrt(fo)) if M>0 else np.arctan2(-dp,-dr)
        def res(z):
            b,s=z
            if not (tlow<=s<=to): return [10,10]
            so,L=ray(M,Xo,b,to,tlow)
            if (so.t[-1]-s)>0: return [5+so.t[-1]-s,5]
            e=em(s)
            if e is None: return [7,7]
            r,p,_=so.sol(s); return cart(r,p)-e[0]
        sol=None
        for guess in ([b0,s0],prev):
            if guess is None: continue
            z,info,ier,msg=fsolve(res,guess,full_output=True,xtol=1e-12)
            if np.max(np.abs(res(z)))<1e-8: sol=z;break
        if sol is None: out.append(np.nan); bad+=1; continue
        prev=sol; b,s=sol; so,L=ray(M,Xo,b,to,tlow); r,p,pr=so.sol(s)
        out.append(kdotu(M,ro,np.cos(b),L,Uo)/kdotu(M,r,pr,L,em(s)[1]))
    return np.array(out),bad
def stats(r): m=np.fft.rfft(r[:-1])/360; return 2*abs(m[1]),2*abs(m[2]),r.max(),r.min()
def norm(d): return d/d.mean()-1

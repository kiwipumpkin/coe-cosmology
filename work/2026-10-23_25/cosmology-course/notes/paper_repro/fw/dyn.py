"""Three frameworks. Units: D (observer distance from CoE at observation) = 1.
CLASSICAL: time unit t_obs=1 (observation time), Newton, static CoE mass mu=G M(<D) t_obs^2/D^3, launch at t_i=0.01 with x=v t_i, v=q*V0.
SR: c=1, time unit D/c. dp/dt=-(K/2) m(x)/x^2, p=gamma v; launch at t_i=1e-4 with proper speed p=q*U. K=2GM(<D)/(D c^2).
GR (LTB): c=1, time unit D/c, shells labelled by areal radius r at t_obs=tau; GM(r)=(K/2) m(r); bang time 0; E(r) from R(r,tau)=r."""
import numpy as np
from core import RG
def _m(x,MT): return np.interp(x,RG,MT)
def logt(ti,tf,NS): s0=np.log(ti); h=(np.log(tf)-s0)/NS; return np.exp(s0+h*np.arange(NS+1)),h
# ---------------- classical ----------------
def cl_integrate(q,V0,mu,MT,NS=400,ti=0.01,tf=1.0,keep=False):
    T,h=logt(ti,tf,NS); x=q*V0*ti; v=q*V0; X=[x];Vv=[v]
    a=lambda x: -mu*_m(x,MT)/np.maximum(x,1e-12)**2
    for n in range(NS):
        t1,t2,t4=T[n],np.sqrt(T[n]*T[n+1]),T[n+1]
        k1=t1*v;l1=t1*a(x); k2=t2*(v+.5*h*l1);l2=t2*a(x+.5*h*k1)
        k3=t2*(v+.5*h*l2);l3=t2*a(x+.5*h*k2); k4=t4*(v+h*l3);l4=t4*a(x+h*k3)
        x=x+h/6*(k1+2*k2+2*k3+k4); v=v+h/6*(l1+2*l2+2*l3+l4)
        if keep: X.append(x);Vv.append(v)
    return (T,np.array(X),np.array(Vv)) if keep else (x,v)
def cl_V0(mu,MT,ti=0.01,NS=400):
    f=lambda V0: cl_integrate(np.array([1.0]),V0,mu,MT,NS=NS,ti=ti)[0][0]-1
    lo,hi=0.5,1.0
    while f(hi)<0: hi*=1.5
    while f(lo)>0: lo/=1.5
    from scipy.optimize import brentq; return brentq(f,lo,hi,xtol=1e-12)
# ---------------- SR ----------------
def sr_integrate(q,U,K,MT,tf,NS=600,ti=1e-4,keep=False):
    T,h=logt(ti,tf,NS); p=q*U; x=p/np.sqrt(1+p*p)*ti; X=[x];P=[p]
    vel=lambda p: p/np.sqrt(1+p*p); f=lambda x: -(K/2)*_m(x,MT)/np.maximum(x,1e-12)**2
    for n in range(NS):
        t1,t2,t4=T[n],np.sqrt(T[n]*T[n+1]),T[n+1]
        k1=t1*vel(p);l1=t1*f(x); k2=t2*vel(p+.5*h*l1);l2=t2*f(x+.5*h*k1)
        k3=t2*vel(p+.5*h*l2);l3=t2*f(x+.5*h*k2); k4=t4*vel(p+h*l3);l4=t4*f(x+h*k3)
        x=x+h/6*(k1+2*k2+2*k3+k4); p=p+h/6*(l1+2*l2+2*l3+l4)
        if keep: X.append(x);P.append(p)
    return (T,np.array(X),np.array(P)) if keep else (x,p)
def sr_tobs(U,K,MT):
    T,X,P=sr_integrate(np.array([1.0]),U,K,MT,tf=1e3,NS=6000,keep=True); x=X[:,0]
    if x.max()<1: return None
    j=np.argmax(x>=1)
    if np.any(P[:j+1,0]<=0): return None
    return float(np.exp(np.interp(0.0,np.log(x[j-1:j+1]),np.log(T[j-1:j+1]))))
# ---------------- GR / LTB ----------------
def _t_of(E,M,r):
    """time for a shell (M,E) to expand from R=0 to R=r (bang time 0), expanding branch; nan if not reachable."""
    E=np.asarray(E,float); out=np.full(np.broadcast(E,M,r).shape,np.nan); E,M,r=np.broadcast_arrays(E,M,r)
    pos=E>1e-300; neg=E<-1e-300; zer=~(pos|neg)
    y=2*E[pos]*r[pos]/np.maximum(M[pos],1e-300); eta=np.arccosh(1+y)
    big=np.where(eta>1e-3,np.sinh(eta)-eta,eta**3/6+eta**5/120)
    out[pos]=np.where(M[pos]>0,M[pos]*big/np.maximum((2*E[pos]),1e-300)**1.5,r[pos]/np.sqrt(2*E[pos]))
    yn=-2*E[neg]*r[neg]/np.maximum(M[neg],1e-300); ok=yn<2
    etn=np.arccos(np.clip(1-yn,-1,1)); sm=np.where(etn>1e-3,etn-np.sin(etn),etn**3/6)
    out[neg]=np.where(ok,M[neg]*sm/(-2*E[neg])**1.5,np.nan)
    out[zer]=(2/3)*np.sqrt(r[zer]**3/(2*np.maximum(M[zer],1e-300)))
    return out
def ltb_E(r,M,tau):
    """solve t(E)=tau per shell by bisection (t decreasing in E)."""
    r=np.asarray(r,float); M=np.asarray(M,float)
    hi=np.full_like(r,1.0); 
    for _ in range(200):
        bad=_t_of(hi,M,r)>tau
        if not bad.any(): break
        hi[bad]*=4
    lo=-M/r*0.999999   # most bound shell reaching r: E>=-M/r
    tl=_t_of(lo,M,r); lo=np.where(np.isfinite(tl)&(tl>=tau),lo,lo)  # if even lo gives t<tau => cannot (would need recollapse); flagged below
    for _ in range(200):
        mid=0.5*(lo+hi); tm=_t_of(mid,M,r); go=np.isfinite(tm)&(tm>tau)
        lo=np.where(go,mid,lo); hi=np.where(go,hi,mid)
    E=0.5*(lo+hi); ok=np.abs(_t_of(E,M,r)/tau-1)<1e-8
    return E,ok
def ltb_R(r,t,Ef,Mf):
    """areal radius R(r,t), Rdot, for shells with E(r),M(r) (callables), bang time 0. vectorised Newton on eta."""
    r=np.asarray(r,float); t=np.asarray(t,float)*np.ones_like(r); E=Ef(r); M=Mf(r)
    R=np.empty_like(r); pos=E>0
    # E>0: R=(M/2E)(cosh eta-1), t=M/(2E)^1.5 (sinh eta-eta)
    A=M[pos]/(2*E[pos]); B=M[pos]/(2*E[pos])**1.5; s=t[pos]/np.maximum(B,1e-300)
    eta=np.where(s>1,np.arcsinh(s)+0.0,np.cbrt(6*s))
    for _ in range(60):
        g=np.sinh(eta)-eta-s; eta=np.maximum(eta-g/(np.cosh(eta)-1+1e-300),1e-12)
    R[pos]=np.where(M[pos]>1e-14*r[pos],A*(np.cosh(eta)-1),np.sqrt(2*E[pos])*t[pos])
    neg=~pos
    if neg.any():
        A=M[neg]/(-2*E[neg]); B=M[neg]/(-2*E[neg])**1.5; s=t[neg]/B; eta=np.cbrt(6*s)
        for _ in range(60):
            g=eta-np.sin(eta)-s; eta=np.clip(eta-g/(1-np.cos(eta)+1e-300),1e-12,np.pi)
        R[neg]=A*(1-np.cos(eta))
    Rd=np.sqrt(np.maximum(2*M/R+2*E,0)); return R,Rd

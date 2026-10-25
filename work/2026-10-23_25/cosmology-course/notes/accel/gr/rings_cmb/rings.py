# Chris's ring dynamics (inferred): each particle moves on its fixed radial line from the CoE (origin), r''=-GM/r^2, r(0)=r0, r'(0)=k r0.
# Matter at t=0 fills the disc |x-(3,0)|<=1. Observer = particle at r0=3 on the axis.
import numpy as np,json,sys
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import brentq
class Flow:
    def __init__(s,k,GM,tmax=6.0):
        s.k,s.GM=k,GM;s.r0=np.linspace(1.9,4.1,221);s.t=np.linspace(0,tmax,601);Rr=[];Vv=[]
        for r0 in s.r0:
            so=solve_ivp(lambda t,y:[y[1],-GM/y[0]**2],(0,tmax),[r0,k*r0],t_eval=s.t,rtol=1e-11,atol=1e-12);Rr.append(so.y[0]);Vv.append(so.y[1])
        s.R=RectBivariateSpline(s.r0,s.t,np.array(Rr),kx=5,ky=5);s.V=RectBivariateSpline(s.r0,s.t,np.array(Vv),kx=5,ky=5)
    def aloc(s,r0,t):  # local 3D scale, axisymmetric about the CoE axis: (dr/dr0 * (r/r0)^2)^(1/3)
        return np.cbrt(s.R(r0,t,dx=1,grid=False)*(s.R(r0,t,grid=False)/r0)**2)
    def r0_of(s,r,t):
        f=lambda q:s.R(q,t,grid=False)-r
        return brentq(f,s.r0[0],s.r0[-1]) if f(s.r0[0])*f(s.r0[-1])<0 else np.nan
def sky(F,c,tobs,arec,N=180,galilean=False):
    ob=F.R(3.0,tobs,grid=False);xo=np.array([ob,0.]);bo=np.array([F.V(3.0,tobs,grid=False)/c,0.]);go=1/np.sqrt(1-bo@bo)
    res=[]
    for th in np.radians(np.arange(N)*360/N):
        # observer-frame sight direction -> lab via aberration
        kp=-np.array([np.cos(th),np.sin(th)]);b=0.0 if galilean else bo[0];kx=(kp[0]+b)/(1+b*kp[0]);ky=kp[1]/(go*(1+b*kp[0]));k=np.array([kx,ky]);k/=np.linalg.norm(k);n=-k
        def g(t):
            x=xo+n*c*(tobs-t);r=np.hypot(*x);q=F.r0_of(r,t)
            if not np.isfinite(q):return np.nan
            p0=q*x/r;inside=np.hypot(p0[0]-3,p0[1])<=1.0
            return (F.aloc(q,t)-arec) if inside else np.nan
        ts=np.linspace(tobs,0,1201);gv=np.array([g(t) for t in ts]);ok=np.isfinite(gv)
        idx=np.where(ok[:-1]&ok[1:]&(np.sign(gv[:-1])!=np.sign(gv[1:])))[0]
        if len(idx)==0:res.append((np.nan,)*5);continue
        i=idx[0];te=brentq(g,ts[i+1],ts[i]);x=xo+n*c*(tobs-te);r=np.hypot(*x);q=F.r0_of(r,te)
        be=F.V(q,te,grid=False)/c*x/r;ge=1/np.sqrt(1-be@be)
        br=be-bo;gr=1/np.sqrt(1-br@br)
        res.append(((1/(gr*(1-br@k))) if galilean else go*(1-bo@k)/(ge*(1-be@k)),te,c*(tobs-te),q,F.aloc(q,te)))
    return np.array(res)

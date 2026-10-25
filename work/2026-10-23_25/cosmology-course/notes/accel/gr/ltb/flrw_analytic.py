# Independent flat-FLRW (matter) check of the no-core ring setup; same initial events/velocities as run_ltb.setup
import sys; sys.path.insert(0,'..'); import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from grcore import relsum, stats, norm
def run(v0,dt,H=0.3,U=0.9999,D=10.0,n=361):
    t1=2/(3*H); a=lambda t:(t/t1)**(2/3); T=t1+dt
    eta=lambda t: 3*t1*(t/t1)**(1/3)          # conformal time, a(t1)=1
    def traj(X0,w,ts):
        g=1/np.sqrt(1-w@w); q=a(ts)*g*w        # comoving momentum, constant
        f=lambda t,X: q/np.sqrt(a(t)**2+q@q)/a(t)
        s=solve_ivp(f,(ts,T),X0,dense_output=True,rtol=1e-12,atol=1e-13)
        vel=lambda t: q/np.sqrt(a(t)**2+q@q)
        return s.sol,vel
    g=1/np.sqrt(1-v0*v0); C0=np.array([D,0.])
    xo,vo=traj(C0+np.array([g*v0*v0,0]),np.array([v0,0.]),t1+g*v0); Xo=xo(T); Vo=vo(T)
    out=[]
    for th in np.radians(np.arange(0,361,360/(n-1))):
        e=np.array([np.cos(th),np.sin(th)]); ts=t1+g*v0*(1+e[0])
        xe,ve=traj(C0+np.array([g*(e[0]+v0*v0),e[1]]),relsum(v0,U,e),ts)
        F=lambda t: np.linalg.norm(Xo-xe(t))-(eta(T)-eta(t))
        te=brentq(F,ts,T); d=Xo-xe(te); nn=d/np.linalg.norm(d); Ve=ve(te)
        go=1/np.sqrt(1-Vo@Vo); ge=1/np.sqrt(1-Ve@Ve)
        out.append(a(te)/a(T)*go*(1-Vo@nn)/(ge*(1-Ve@nn)))
    return np.array(out)
if __name__=='__main__':
    for v0,dt in [(0,1.4),(0.1,1.6),(0.3,1.9),(0.5,2.4),(0.9,10.5)]:
        d,q,_,mn=stats(norm(run(v0,dt))); print(v0,"%.4f %.4f %.4f"%(d,q,mn),flush=True)

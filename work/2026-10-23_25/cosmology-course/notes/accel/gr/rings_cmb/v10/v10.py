# v10: Chris's kinematics. Shells r(t,w): dr/dt=sqrt(2M/r+2E), E=w|w|/2, common start t=0 (Newtonian time, speeds may exceed c).
# Light: dx/dt = v(r,t) rhat + c n  (c relative to local matter); direction n turned by the velocity shear;
# redshift: d ln(1+z) = n.(grad v).n dt  = [v' cos^2a + (v/r) sin^2a] dt  (accumulated Doppler between neighbouring matter).
import sys,numpy as np,v5 as V,pickle
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline
uD=10**0.021927885054513485
def field(M):
    wo,t0=V.setup(M)
    U=np.r_[np.linspace(np.arcsinh(-0.9 if M>0 else 1e-4),np.arcsinh(2),600),np.linspace(np.arcsinh(2),np.arcsinh(60),300)[1:]]
    T=np.geomspace(1e-5,t0,700)
    G=np.log(np.maximum([[V.Rfun(t,np.sinh(u),M) for u in U] for t in T],1e-14))
    return wo,t0,U,T,RectBivariateSpline(np.log(T),U,G,kx=5,ky=5)
def rays(M,Rrec,psis,F=None):
    wo,t0,U,T,S=F or field(M)
    def q(t,u):
        lt=np.log(t);L=S(lt,u)[0,0];R=np.exp(L);Lu=S(lt,u,dy=1)[0,0];Lt=S(lt,u,dx=1)[0,0]/t;Ltu=S(lt,u,dx=1,dy=1)[0,0]/t
        Ru=R*Lu;v=R*Lt;vu=R*(Ltu+Lt*Lu);return R,Ru,v,vu/Ru
    def rhs(t,y):
        u,th,a,L=y;R,Ru,v,vp=q(t,u);c,s=np.cos(a),np.sin(a)
        return [c/Ru, s/R, s*c*(vp-v/R)-s/R, vp*c*c+(v/R)*s*s]
    hit=lambda t,y:q(t,y[0])[0]-Rrec; hit.terminal=True
    edge=lambda t,y:min(y[0]-U[0]-1e-3,U[-1]-1e-3-y[0]); edge.terminal=True
    uo=np.arcsinh(wo);out=[];paths=[]
    for ps in psis:
        y0=[uo,0.,np.radians(ps)+np.pi,0.]
        s=solve_ivp(rhs,(t0,T[0]*1.01),y0,events=[hit,edge],rtol=1e-9,atol=1e-12,dense_output=True)
        tt=np.linspace(t0,s.t[-1],200);Y=s.sol(tt);paths.append((ps,tt,np.array([q(a,b)[0] for a,b in zip(tt,Y[0])]),Y[1]))
        if s.t_events[0].size:
            te=s.t_events[0][0];y=s.y_events[0][0];out.append((ps,np.exp(-y[3]) and np.exp(-y[3]),te,t0-te,np.sinh(y[0]),y[1],0))
        else:
            y=s.y[:,-1];code=1 if s.t_events[1].size else 2;out.append((ps,np.nan,s.t[-1],np.nan,np.sinh(y[0]),y[1],code))
    return np.array(out),paths,wo,t0
if __name__=='__main__':
    M=float(sys.argv[1]);F=field(M);res={}
    for R in [float(x) for x in sys.argv[2].split(',')]:
        o,p,wo,t0=rays(M,R,np.arange(0,181,1.0),F);res[R]=(o,p);print(M,R,'hits',int(np.isfinite(o[:,1]).sum()),'codes',np.unique(o[:,6],return_counts=True),flush=True)
    pickle.dump(dict(res=res,wo=wo,t0=t0),open(f'v10_M{M:.5f}_{sys.argv[3]}.pkl','wb'))

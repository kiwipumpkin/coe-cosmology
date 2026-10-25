# v11: comoving coords (label w, phi), metric ds^2=-dt^2+(dr/dw)^2 dw^2+r^2 dphi^2 ; r(t,w) from shell solution dr/dt=sqrt(2M/r+w|w|), common start.
import sys,numpy as np,base as V,pickle
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline
def field(M,Rmin):
    wo,t0=V.setup(M);wmin=-np.sqrt(2*M/Rmin)*0.98 if M>0 else 1e-3
    T=np.geomspace(2e-3,t0,500);W=np.r_[np.linspace(wmin,1.5,500),np.geomspace(1.5,40,200)[1:]]
    S=RectBivariateSpline(T,W,np.array([[V.Rfun(t,w,M) for w in W] for t in T]),kx=5,ky=5);return wo,t0,T,W,S
def trace(F,Rrec,psis,tobs=None,keep=False):
    wo,t0,T,W,S=F;tobs=tobs or t0
    def AB(t,w):
        R=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0]
        return R,Rp*Rp,2*Rp*Rtp,2*Rp*Rpp,R*R,2*R*Rt,2*R*Rp
    def rhs(l,y):
        t,w,ph,pt,pw,pp=y;R,A,At,Aw,B,Bt,Bw=AB(t,w)
        return [-pt,pw/A,pp/B,0.5*At/A**2*pw**2+0.5*Bt/B**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B**2*pp**2,0]
    hit=lambda l,y:AB(y[0],y[1])[0]-Rrec;hit.terminal=True
    edge=lambda l,y:min(y[0]-T[0]*1.01,y[1]-W[0]-1e-3,W[-1]-1e-2-y[1]);edge.terminal=True
    A0=AB(tobs,wo)[1];out=[];paths=[]
    for ps in psis:
        pr=np.radians(ps);y0=[tobs,wo,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sqrt(AB(tobs,wo)[4])*np.sin(pr)]
        s=solve_ivp(rhs,(0,-200),y0,events=[hit,edge],rtol=1e-9,atol=1e-11,dense_output=keep,max_step=(0.01 if keep else np.inf))
        if keep: paths.append((ps,s.y[0],s.y[1],s.y[2]))
        if s.t_events[0].size:
            y=s.y_events[0][0];out.append((ps,-y[3],y[0],tobs-y[0],y[1],y[2],0))
        else: y=s.y[:,-1];out.append((ps,np.nan,y[0],np.nan,y[1],y[2],1))
    return np.array(out),paths
if __name__=='__main__':
    M=float(sys.argv[1]);Rs=[float(x) for x in sys.argv[2].split(',')];F=field(M,min(Rs));res={}
    for R in Rs:
        o,_=trace(F,R,np.arange(0,181,2.0));res[R]=o;print(M,R,int(np.isfinite(o[:,1]).sum()),flush=True)
    pickle.dump(dict(res=res,t0=F[1],wo=F[0]),open(f'v11_M{M:.5f}_{sys.argv[3]}.pkl','wb'))

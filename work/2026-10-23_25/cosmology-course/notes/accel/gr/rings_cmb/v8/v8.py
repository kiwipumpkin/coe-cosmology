# v8: emission event set purely by CoE distance: element emits when its proper (common-bang) time tau = t_e(R) = tau_s * R^(-q)  (R in units of D; outer first).
import sys,numpy as np,v5 as V
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline
M=float(sys.argv[1]); fr=float(sys.argv[2]); q=float(sys.argv[3]); a,b=float(sys.argv[4]),float(sys.argv[5])
wo,t0=V.setup(M); ts=fr*V.setup(0)[1]; tmin=0.2*ts*min(1,0.05**(-q)) if False else 0.15*ts; wmax=max(60,8*V.setup(0)[1]/ts)
wmin=-np.sqrt(2*M/0.02) if M>0 else 1e-3
T=np.geomspace(tmin,t0,600); W=np.r_[np.linspace(wmin,1.5,400),np.geomspace(1.5,wmax,400)[1:]]
S=RectBivariateSpline(T,W,np.array([[V.Rfun(t,w,M) for w in W] for t in T]),kx=5,ky=5)
E=lambda w:np.sign(w)*w*w/2
def AB(t,w):
    R=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0];qq=1+2*E(w)
    return R,Rp**2/qq,2*Rp*Rtp/qq,2*Rp*Rpp/qq-Rp**2*2*abs(w)/qq**2,R*R,2*R*Rt,2*R*Rp
def rhs(l,y):
    t,w,ph,pt,pw,pp=y;R,A,At,Aw,B,Bt,Bw=AB(t,w)
    return [-pt,pw/A,pp/B,0.5*At/A**2*pw**2+0.5*Bt/B**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B**2*pp**2,0]
def hit(l,y): R=max(AB(y[0],y[1])[0],1e-9); return np.log(y[0])-np.log(ts)+q*np.log(R)
hit.terminal=True
edge=lambda l,y:min(y[0]-T[0]*1.001,y[1]-W[0]-1e-3,W[-1]*0.999-y[1]);edge.terminal=True
A0=AB(t0,wo)[1];out=[]
for ps in np.arange(a,b+0.5,2.0):
    pr=np.radians(ps);y0=[t0,wo,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sin(pr)]
    s=solve_ivp(rhs,(0,-500),y0,events=[hit,edge],rtol=1e-9,atol=1e-12)
    if s.t_events[0].size: y=s.y_events[0][0];out.append((ps,-y[3],y[0],t0-y[0],y[1],AB(y[0],y[1])[0]))
    else: out.append((ps,)+(np.nan,)*5)
np.save(f'p_M{M:.5f}_f{fr}_q{q}_{int(a)}.npy',np.array(out));print(wo,t0,ts,np.isfinite(np.array(out)[:,1]).sum(),flush=True)

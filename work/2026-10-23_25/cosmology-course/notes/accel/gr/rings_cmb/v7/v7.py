# v6: recombination when each element's own plasma density (T ∝ n^{1/3}, adiabatic) falls to n_rec.
# n ∝ w^2/(R^2 R') (label w = asymptotic-speed parameter, E=w|w|/2); this is exactly homogeneous Milne (n=1/tau^3) when M=0.
import sys,numpy as np,v5 as V
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline
def build(M,t0,wmin,wmax=60):
    T=np.geomspace(5e-3,t0,500); W=np.r_[np.linspace(wmin,1.5,400),np.geomspace(1.5,wmax,300)[1:]]
    G=np.array([[V.Rfun(t,w,M) for w in W] for t in T]); return RectBivariateSpline(T,W,G,kx=5,ky=5),T,W
def trace(M,taurec,psis,p=0.0):
    wo,t0=V.setup(M); wmin=-np.sqrt(2*M/0.05) if M>0 else 1e-3
    S,T,W=build(M,t0,wmin); E=lambda w:np.sign(w)*w*w/2; wo0=V.setup(0)[0]
    def AB(t,w):
        R=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0];q=1+2*E(w)
        A=Rp**2/q; At=2*Rp*Rtp/q; Aw=2*Rp*Rpp/q-Rp**2*2*abs(w)/q**2
        return R,Rp,A,At,Aw,R*R,2*R*Rt,2*R*Rp
    def rhs(l,y):
        t,w,ph,pt,pw,pp=y; R,Rp,A,At,Aw,B,Bt,Bw=AB(t,w)
        return [-pt,pw/A,pp/B,0.5*At/A**2*pw**2+0.5*Bt/B**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B**2*pp**2,0]
    def hit(l,y):
        R,Rp,*_=AB(y[0],y[1]); return np.log(y[1]**2/(R*R*abs(Rp))+1e-300)+3*np.log(taurec)-3*p*np.log(abs(y[1])/wo0+1e-300)
    hit.terminal=True
    edge=lambda l,y:min(y[0]-T[0]-1e-4,y[1]-W[0]-1e-3,W[-1]-y[1]-1e-2); edge.terminal=True
    _,_,A0,*_=AB(t0,wo); n0=wo**2/(1*abs(AB(t0,wo)[1])); out=[]
    for ps in psis:
        pr=np.radians(ps); y0=[t0,wo,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sin(pr)]
        s=solve_ivp(rhs,(0,-200),y0,events=[hit,edge],rtol=1e-9,atol=1e-11)
        if s.t_events[0].size: y=s.y_events[0][0]; out.append((-y[3],y[0],t0-y[0],y[1]))
        else: out.append((np.nan,)*4)
    return np.array(out),wo,t0,n0
if __name__=='__main__':
    M=float(sys.argv[1]); fr=float(sys.argv[2]); a,b=float(sys.argv[3]),float(sys.argv[4]); P=float(sys.argv[5])
    wo,t0=V.setup(M); tr=fr*V.setup(0)[1]       # tau_rec = fr * (no-mass age); same n_rec for both cases
    r,wo,t0,n0=trace(M,tr,np.arange(a,b+0.5,1.0),P)
    np.save(f'p_M{M:.5f}_f{fr}_p{P}_{int(a)}.npy',np.c_[np.arange(a,b+0.5,1.0),r]); print(wo,t0,n0,tr,np.isfinite(r[:,0]).sum(),flush=True)

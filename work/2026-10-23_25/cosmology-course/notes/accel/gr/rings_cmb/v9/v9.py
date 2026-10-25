# v9: computed cooling. Initial state: uniform proper density and uniform temperature T_i at early proper time tau_i (same both cases).
# n/n_i = R_i^2 dR_i/du / (R^2 dR/du) (shell solution incl. mass), T = T_i (n/n_i)^(1/3)  (radiation-dominated, gamma=4/3).
# Emits when T = 3000 K. Normalisation fixed so the no-mass case has 1+z = Z0 exactly.
import sys,numpy as np,v5 as V
from scipy.integrate import solve_ivp
from scipy.interpolate import RectBivariateSpline,CubicSpline
M=float(sys.argv[1]); Z0=float(sys.argv[2]); ti_frac=float(sys.argv[3]); a,b=float(sys.argv[4]),float(sys.argv[5])
wo,t0=V.setup(M); t0n=V.setup(0)[1]; te0=t0n/Z0; ti=ti_frac*te0; lnX=3*np.log(te0/ti)   # (T_i/3000)^3
wmin=-0.6 if M>0 else 1e-4; wmax=40*Z0
U=np.r_[np.linspace(np.arcsinh(wmin),np.arcsinh(2.0),500),np.linspace(np.arcsinh(2.0),np.arcsinh(wmax),500)[1:]]; Wg=np.sinh(U)
T=np.geomspace(0.25*te0,t0,700)
G=np.log(np.maximum([[V.Rfun(t,w,M) for w in Wg] for t in T],1e-12))
S=RectBivariateSpline(np.log(T),U,G,kx=5,ky=5)
import mpmath as mp; mp.mp.dps=50
def Rmp(tau,w):
    tau=mp.mpf(tau);w=mp.mpf(w);E=w*abs(w)/2;Mm=mp.mpf(M)
    if M==0: return abs(w)*tau
    a=Mm/(2*abs(E));b=a/mp.sqrt(2*abs(E));x=tau/b
    if E>0: eta=mp.findroot(lambda e:mp.sinh(e)-e-x,mp.cbrt(6*x) if x<1 else mp.log(2*x+1)+1);return a*(mp.cosh(eta)-1)
    eta=mp.findroot(lambda e:e-mp.sin(e)-x,mp.cbrt(6*x));return a*(1-mp.cos(eta))
def dlnRi(w):
    if M==0: return 1/w
    h=mp.mpf(1e-12)*max(abs(w),1e-3); return float((mp.log(Rmp(ti,w+h))-mp.log(Rmp(ti,w-h)))/(2*h))
LRi=np.array([float(mp.log(Rmp(ti,w))) for w in Wg]); dLi=np.array([dlnRi(w) for w in Wg])*np.cosh(U)
Si=CubicSpline(U,LRi); Sd=CubicSpline(U,np.log(np.abs(dLi)))
def init(u): return 3*Si(u)+Sd(u)      # ln(R_i^2 dR_i/du) up to const (R_i' = R_i * dlnR_i/du)
E=lambda w:np.sign(w)*w*w/2
def AB(t,u):
    lt=np.log(t);L=S(lt,u)[0,0];Lu=S(lt,u,dy=1)[0,0];Lt=S(lt,u,dx=1)[0,0]/t;Luu=S(lt,u,dy=2)[0,0];Ltu=S(lt,u,dx=1,dy=1)[0,0]/t
    R=np.exp(L);Ru=R*Lu;Ruu=R*(Luu+Lu*Lu);Rt=R*Lt;Rtu=R*(Ltu+Lt*Lu);w=np.sinh(u);q=1+2*E(w);dq=2*abs(w)*np.cosh(u)
    return R,Ru,Ru**2/q,2*Ru*Rtu/q,2*Ru*Ruu/q-Ru**2*dq/q**2,R*R,2*R*Rt,2*R*Ru,L,Lu
def rhs(l,y):
    t,u,ph,pt,pu,pp=y;R,Ru,A,At,Au,B,Bt,Bu,*_=AB(t,u)
    return [-pt,pu/A,pp/B,0.5*At/A**2*pu**2+0.5*Bt/B**2*pp**2,0.5*Au/A**2*pu**2+0.5*Bu/B**2*pp**2,0]
def lnn(t,u):  # ln(n/n_i)
    *_,L,Lu=AB(t,u); return init(u)-(3*L+np.log(abs(Lu)))
hit=lambda l,y:lnn(y[0],y[1])+lnX; hit.terminal=True
edge=lambda l,y:min(np.log(y[0]/T[0]),y[1]-U[0]-1e-3,U[-1]-1e-3-y[1]);edge.terminal=True
uo=np.arcsinh(wo);A0=AB(t0,uo)[2];out=[]
for ps in np.arange(a,b+0.5,2.0):
    pr=np.radians(ps);y0=[t0,uo,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sin(pr)]
    l0=0;nref=0
    while True:
        s=solve_ivp(rhs,(l0,l0-1e4),y0,events=[hit,edge],rtol=1e-10,atol=1e-13)
        if M==0 and s.t_events[1].size and s.y_events[1][0][1]<U[0]+0.01 and nref<3:   # pass through regular centre
            y0=list(s.y_events[1][0]);y0[1]=U[0]+0.0011;y0[4]=-y0[4];y0[2]+=np.pi;l0=s.t_events[1][0];nref+=1;continue
        break
    if s.t_events[0].size:
        y=s.y_events[0][0];w=np.sinh(y[1]);fb=float(w<0 and y[0]>np.pi/2*(2*M/w**2)**1.5/np.sqrt(2*M))   # past turnaround?
        out.append((ps,-y[3],y[0],t0-y[0],w,AB(y[0],y[1])[0],fb,0))
    else:
        y=s.y[:,-1];out.append((ps,np.nan,y[0],np.nan,np.sinh(y[1]),np.nan,np.nan,1))
np.save(f'p_M{M:.5f}_Z{int(Z0)}_ti{ti_frac}_{int(a)}.npy',np.array(out));print(wo,t0,te0,ti,np.isfinite(np.array(out)[:,1]).sum(),flush=True)

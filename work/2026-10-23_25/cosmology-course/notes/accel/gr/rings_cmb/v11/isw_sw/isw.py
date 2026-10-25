# Order-of-magnitude Newtonian SW + ISW estimate along the model's own ray-traced LTB geodesics.
# Units: G=c=1, D(us today)=1. Phi(t,x) = -M(<x)/R(t,x) - sum_{x'>x} dM(x')/R(t,x')  (Newtonian, spherical shells, no shell crossing).
import sys,os,json,numpy as np
fam,p,frac,out=sys.argv[1],float(sys.argv[2]),float(sys.argv[3]),sys.argv[4]
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/';sys.path.insert(0,B+'dipole_amp');os.environ.setdefault('TMIN','1e-8');os.chdir(B+'dipole_amp')
import dist2 as D
from scipy.optimize import brentq
from scipy.interpolate import interp1d,RectBivariateSpline
from scipy.integrate import solve_ivp,cumulative_trapezoid
SLOPE=0.04103732750534994;DR=42/14.30048287260589;Rb=1.5*(1+1/DR)
base=lambda x:1/(1+np.abs(x)**p)
ls=lambda L:-(np.log(base(1.001/L))-np.log(base(0.999/L)))/0.002-SLOPE
g=np.geomspace(0.05,1e5,600);v=np.array([ls(L) for L in g]);i=np.where(np.sign(v[:-1])!=np.sign(v[1:]))[0][-1]
L=brentq(ls,g[i],g[i+1]);shape=lambda x:base(x/L)
X=np.concatenate([[0.0],np.geomspace(1e-5,60,6000)]);sv=shape(np.maximum(X,1e-12))*(X<=Rb)
m=cumulative_trapezoid(4*np.pi*X**2*sv,X,initial=0.0);mt=(X,m/np.interp(1.0,X,m))
F=D.build(frac*0.5*D.u**2,mt);S,t0,XL,Mv,T=F['S'],F['t0'],F['X'],F['M'],F['T']
# potential grid
G=F['G'];dM=np.gradient(Mv,XL)
PHI=np.empty_like(G)
for k in range(len(T)):
    Rk=G[k];outer=np.r_[cumulative_trapezoid((dM/Rk)[::-1],XL[::-1],initial=0)[::-1]]*-1  # int_x^Xmax dM/R
    PHI[k]=-Mv/Rk-outer
lT=np.log(T);SP=RectBivariateSpline(lT,XL,PHI,kx=3,ky=3)
Phi=lambda t,x:SP(np.log(t),x)[0,0];dPhidt=lambda t,x:SP(np.log(t),x,dx=1)[0,0]/t
def rho(t,x):
    R=S(t,x)[0,0];Rp=S(t,x,dy=1)[0,0];return shape(x)*x*x/(R*R*Rp)
tus=brentq(lambda t:S(t,1.0)[0,0]-1/1100.,T[0]*1.001,t0);rc=rho(tus,1.0)
XX=np.linspace(0.6,1.6,201);te=np.array([brentq(lambda t:np.log(rho(t,x)/rc),T[0]*1.0001,t0,xtol=1e-16,rtol=1e-14) for x in XX]);tf=interp1d(XX,te,fill_value='extrapolate')
def AB(t,w):
    Rr=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0]
    return Rr,Rp*Rp,2*Rp*Rtp,2*Rp*Rpp,Rr*Rr,2*Rr*Rt,2*Rr*Rp
def rhs(l,y):
    t,w,ph,pt,pw,pp,I=y;Rr,A,At,Aw,B_,Bt,Bw=AB(t,w)
    return [-pt,pw/A,pp/B_,0.5*At/A**2*pw**2+0.5*Bt/B_**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B_**2*pp**2,0,2*dPhidt(t,w)*(-pt)]  # dI = 2 dPhi/dt dt
hit=lambda l,y:y[0]-float(tf(y[1]));hit.terminal=True
A0=AB(t0,1.0)[1];B0=AB(t0,1.0)[4];o=[]
for ps in np.arange(0,180.001,5):
    pr=np.radians(ps);y0=[t0,1.0,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sqrt(B0)*np.sin(pr),0.]
    s=solve_ivp(rhs,(0,-200),y0,events=[hit],rtol=1e-10,atol=1e-13);y=s.y_events[0][0]
    # I accumulated backward in t: ISW = 2∫_{te}^{t0} dPhi/dt dt = -I_final (since dt<0 going back)
    o.append((ps,-y[3],y[0],y[1],Phi(y[0],y[1]),-y[6]))
o=np.array(o);np.savez(out,o=o,frac=frac,p=p,L=L,Phi_obs=Phi(t0,1.0),t0=t0)
print(out,'1+z',o[0,1],o[-1,1],'Phi_e',o[:,4].min(),o[:,4].max(),'ISW',o[:,5].min(),o[:,5].max(),'Phi_obs',Phi(t0,1.0),flush=True)

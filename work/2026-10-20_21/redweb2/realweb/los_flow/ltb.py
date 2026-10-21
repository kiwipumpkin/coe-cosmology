# Red-universe LTB (famP 'poly' build, frac=20000, u=4, D=42 Gpc): expansion law + photon paths in Lagrangian coords + growth factors.
import sys,os,numpy as np
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/';sys.path.insert(0,B+'dipole_amp');os.environ.setdefault('TMIN','1e-8');os.chdir(B+'dipole_amp')
import dist2 as D
from scipy.integrate import solve_ivp,cumulative_trapezoid
from scipy.optimize import brentq
from scipy.interpolate import interp1d
DR=42/14.30048287260589;Rb=1.5*(1+1/DR);CC=[-0.041,-0.035,0.1353,0.0];frac=20000.
shape=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)))
X=np.concatenate([[0.0],np.geomspace(1e-5,60,6000)]);sv=shape(np.maximum(X,1e-12))*(X<=Rb)
m=cumulative_trapezoid(4*np.pi*X**2*sv,X,initial=0.0);mt=(X,m/np.interp(1.0,X,m))
F=D.build(frac*0.5*D.u**2,mt);S,t0=F['S'],F['t0']
def rho(t,x):R=S(t,x)[0,0];Rp=S(t,x,dy=1)[0,0];return shape(x)*x*x/(R*R*Rp)
tus=brentq(lambda t:S(t,1.0)[0,0]-1/1100.,F['T'][0]*1.001,t0);rc=rho(tus,1.0)
XX=np.linspace(0.98,1.02,81);te=np.array([brentq(lambda t:np.log(rho(t,x)/rc),F['T'][0]*1.0001,t0,xtol=1e-16,rtol=1e-14) for x in XX]);tf=interp1d(XX,te,fill_value='extrapolate')
def AB(t,w):
    Rr=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0]
    return Rr,Rp*Rp,2*Rp*Rtp,2*Rp*Rpp,Rr*Rr,2*Rr*Rt,2*Rr*Rp
def rhs(l,y):
    t,w,ph,pt,pw,pp=y;Rr,A,At,Aw,B_,Bt,Bw=AB(t,w)
    return [-pt,pw/A,pp/B_,0.5*At/A**2*pw**2+0.5*Bt/B_**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B_**2*pp**2,0]
hit=lambda l,y:y[0]-float(tf(y[1]));hit.terminal=True
A0=AB(t0,1.0)[1];B0=AB(t0,1.0)[4];PS=np.arange(0,180.001,5);paths=[]
for ps in PS:
    pr=np.radians(ps);y0=[t0,1.0,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sqrt(B0)*np.sin(pr)]
    s=solve_ivp(rhs,(0,-200),y0,events=[hit],rtol=1e-10,atol=1e-13,dense_output=True);le=s.t_events[0][0]
    # sample uniformly in log(t0 - ... ) : use 400 points in lambda, plus log-spaced near source (time compresses)
    L=np.r_[np.linspace(0,le,300)];Y=s.sol(L);paths.append(np.stack([L,Y[0],Y[1],Y[2]],1))
    print(ps,'x_src',Y[1,-1],'phi_src',Y[2,-1],'t_src',Y[0,-1],flush=True)
# expansion and growth at our shell x=1 (local, isotropic approx; anisotropy q=(R'/R'0)/(R/R0) reported)
T=np.geomspace(F['T'][0]*1.01,t0,3000);a_t=np.array([S(t,1.)[0,0] for t in T])/S(t0,1.)[0,0]
a_r=np.array([S(t,1.,dy=1)[0,0] for t in T])/S(t0,1.,dy=1)[0,0];q=a_r/a_t
Hh=np.array([S(t,1.,dx=1)[0,0]/S(t,1.)[0,0] for t in T])
Mp=F['Mtot']*np.gradient(np.interp(np.array([0.999,1.0,1.001]),*mt),0.001)[1]
fourpiGrho=np.array([Mp/(S(t,1.)[0,0]**2*S(t,1.,dy=1)[0,0]) for t in T])   # 4 pi G rho = M'(x)/(R^2 R')
# growth: d''+2H d'=4piGrho d, in s=ln t ; start growing mode d ~ t^(2/3)
fr=interp1d(np.log(T),fourpiGrho);fH=interp1d(np.log(T),Hh)
def g(s,y):t=np.exp(s);d,dd=y;return [t*dd,t*(-2*fH(s)*dd+fr(s)*d)]
s0=np.log(T[0]);sol=solve_ivp(g,(s0,np.log(T[-1])),[1.0,(2/3)/T[0]],t_eval=np.log(T),rtol=1e-9,atol=1e-14,max_step=0.01)
Dg=sol.y[0]/sol.y[0][-1];Agrow=Dg/a_t;Afix=1/a_t
np.savez('/workspace/redweb2/realweb/los_flow/ltb.npz',T=T,t0=t0,a_t=a_t,a_r=a_r,q=q,H=Hh,Agrow=Agrow/Agrow[-1],Afix=Afix,PS=PS,paths=np.array(paths),Omega=2*fourpiGrho/3/Hh**2)
print('t0',t0,'t_rec~',T[0],'q range',q.min(),q.max(),'Agrow at rec',Agrow[0]/Agrow[-1],'min/max',Agrow.min()/Agrow[-1],Agrow.max()/Agrow[-1],'Omega now',(2*fourpiGrho/3/Hh**2)[-1])

# Shrinking-radius (critical density) CMB trace, 1+z~1100, for a density family re-locked to the CatWISE log-slope -0.041037 at our distance.
# Same shape is used for the mass distribution (enclosed-mass table, truncated at Rb as before) and the plasma density.
import sys,os,json,numpy as np
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/';sys.path.insert(0,B+'dipole_amp');os.environ.setdefault('TMIN','1e-8');os.chdir(B+'dipole_amp')
import dist2 as D
from scipy.optimize import brentq
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp,cumulative_trapezoid
fam,p,frac,out=sys.argv[1],float(sys.argv[2]),float(sys.argv[3]),os.path.abspath(sys.argv[4]) if sys.argv[4].startswith('/') else sys.argv[4]
SLOPE=0.04103732750534994;DR=float(os.environ.get("DRV",42/14.30048287260589));Rb=1.5*(1+1/DR)
def base(kind,x,p):
    x=np.abs(x)
    return {'soft':lambda:1/(1+x**p),'exp':lambda:np.exp(-x),'eta':lambda:np.exp(-x**p)}[kind]()
if fam=='curv':  # slope locked at x=1, far/core curvature free: ln rho = -SLOPE ln x - p (x-1)^2
    shape=lambda x:np.abs(x)**-SLOPE*np.exp(-p*(np.abs(x)-1)**2);L=None
elif fam=='cub':  # slope locked: ln rho = -SLOPE ln x - p (x-1)^2 - q (x-1)^3, q from env Q
    qq=float(os.environ.get('Q','0'));shape=lambda x:np.abs(x)**-SLOPE*np.exp(-p*(np.abs(x)-1)**2-qq*(np.minimum(np.abs(x),Rb)-1)**3);L=None
elif fam=='slc':  # FREE near-slope SL (env SL) + curvature p: ln rho = -SL ln x - p (x-1)^2
    SLv=float(os.environ.get('SL','0'));shape=lambda x:np.abs(x)**-SLv*np.exp(-p*(np.abs(x)-1)**2);L=None
elif fam=='outer':  # slope locked at x=1, mass rises beyond observer: ln rho = -SLOPE ln x + p*max(x-1,0)^2
    shape=lambda x:np.abs(x)**-SLOPE*np.exp(p*np.maximum(np.abs(x)-1,0)**2);L=None
elif fam=='slope':  # unlocked power law rho ~ x^-p (p=0.041 is the lock; p<=0 is outer-weighted)
    shape=lambda x:np.abs(x)**-p;L=None
elif fam=='core':  # slope locked at x=1 (core term ~e^-10 there), add central overdensity p
    shape=lambda x:np.abs(x)**-SLOPE*(1+p*np.exp(-np.abs(x)/0.05));L=None
elif fam=='pow': shape=lambda x:np.abs(x)**-SLOPE if True else 0;L=None
else:
    ls=lambda L:-(np.log(base(fam,1.001/L,p))-np.log(base(fam,0.999/L,p)))/0.002-SLOPE
    g=np.geomspace(0.05,1e5,600);v=np.array([ls(L) for L in g]);i=np.where(np.sign(v[:-1])!=np.sign(v[1:]))[0][-1]
    L=brentq(ls,g[i],g[i+1]);shape=lambda x:base(fam,x/L,p)
X=np.concatenate([[0.0],np.geomspace(1e-5,60,6000)]);sv=np.where(X>0,shape(np.maximum(X,1e-12)),shape(1e-12))*(X<=Rb)
m=cumulative_trapezoid(4*np.pi*X**2*sv,X,initial=0.0);mt=(X,m/np.interp(1.0,X,m))
F=D.build(frac*0.5*D.u**2,mt);S,t0=F['S'],F['t0']
R1=S(t0,1.)[0,0];Rp=S(t0,1.,dy=1)[0,0];Rt=S(t0,1.,dx=1)[0,0];Rtp=S(t0,1.,dx=1,dy=1)[0,0];ratio=(Rtp/Rp)/(Rt/R1)
def rho(t,x):
    R=S(t,x)[0,0];Rp=S(t,x,dy=1)[0,0];return shape(x)*x*x/(R*R*Rp)
tus=brentq(lambda t:S(t,1.0)[0,0]-1/1100.,F['T'][0]*1.001,t0);rc=rho(tus,1.0)
XX=np.linspace(0.6,1.6,201);te=np.array([brentq(lambda t:np.log(rho(t,x)/rc),F['T'][0]*1.0001,t0,xtol=1e-16,rtol=1e-14) for x in XX])
tf=interp1d(XX,te,fill_value='extrapolate')
def AB(t,w):
    Rr=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0]
    return Rr,Rp*Rp,2*Rp*Rtp,2*Rp*Rpp,Rr*Rr,2*Rr*Rt,2*Rr*Rp
def rhs(l,y):
    t,w,ph,pt,pw,pp=y;Rr,A,At,Aw,B_,Bt,Bw=AB(t,w)
    return [-pt,pw/A,pp/B_,0.5*At/A**2*pw**2+0.5*Bt/B_**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B_**2*pp**2,0]
hit=lambda l,y:y[0]-float(tf(y[1]));hit.terminal=True
A0=AB(t0,1.0)[1];B0=AB(t0,1.0)[4];o=[];P=[]
for ps in [0,20,45,70,90,110,135,160,180]:
    pr=np.radians(ps);y0=[t0,1.0,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sqrt(B0)*np.sin(pr)]
    s=solve_ivp(rhs,(0,-200),y0,events=[hit],rtol=1e-10,atol=1e-13,dense_output=True);y=s.y_events[0][0]
    lg=np.linspace(0,s.t_events[0][0],2001);Y=s.sol(lg);tt=Y[0];ww=Y[1];Q=np.array([AB(a,b) for a,b in zip(tt,ww)]);E=-Y[3]
    c2=(Y[4]/np.sqrt(Q[:,1])/E)**2;s2=(Y[5]/np.sqrt(Q[:,4])/E)**2;Hr=Q[:,2]/(2*Q[:,1]);Hp=Q[:,5]/(2*Q[:,4]);P.append(dict(psi=ps,w=ww,ph=Y[2],R=np.array([S(a,b)[0,0] for a,b in zip(tt[::20],ww[::20])]),t=tt,c2=c2,s2=s2,Hr=Hr,Hp=Hp,lnE=np.log(E)))
    o.append((ps,-y[3],y[0],t0-y[0],y[1],y[2],S(y[0],y[1])[0,0]))
import pickle;pickle.dump(P,open(out+'.xy.pkl','wb'));np.savez(out,o=np.array(o),fam=fam,p=p,frac=frac,L=L if L else np.nan,ratio=ratio,q=float(os.environ.get('Q','0')),sl=float(os.environ.get('SL','nan')),te=te,XX=XX)
print(out,'L/D',L,'ratio',ratio,'1+z',o[0][1],o[-1][1],flush=True)

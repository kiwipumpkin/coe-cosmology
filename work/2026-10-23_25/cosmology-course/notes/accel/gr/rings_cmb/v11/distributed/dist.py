# Distributed-mass (soft core+tail) version of the v11 comoving light tracing at 4c, D = 42 Gpc, R_rec = 0.5 D.
# Units: D = 1 (our distance today), c = 1, G = 1. Our speed today u = 4.
# Mass: M(<x) = Mtot * m(x), m(x) = enclosed-mass fraction of rho ∝ 1/(1+(x/LD)^alpha), truncated at Rb (as in prof_core), m(1) = 1.
# Time dependence: the law is defined on TODAY's radii. Each shell (label x = its radius today) keeps its enclosed mass M(<x) for all time
# (Newtonian shells, Gauss's law), valid as long as shells never cross (checked). Common bang: every shell starts at r=0 at t=0.
import sys,json,numpy as np
sys.path.insert(0,'..');import base as V
from scipy.integrate import quad,cumulative_trapezoid
from scipy.optimize import brentq
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import solve_ivp
R=0.5;u=4.0
def mfrac_table(LD,alpha,Rb,kind='soft'):
    X=np.concatenate([[0.0],np.geomspace(1e-5,60,6000)])
    if kind=='point': return X,np.ones_like(X)
    s=(1+np.abs(X/LD)**alpha)**-1.0*(X<=Rb)
    m=cumulative_trapezoid(4*np.pi*X**2*s,X,initial=0.0);return X,m/np.interp(1.0,X,m)
def age_out(r,M,E):
    if M==0: return r/np.sqrt(2*E)
    return quad(lambda y:1/np.sqrt(max(2*M/y+2*E,1e-300)),0,r,limit=200)[0]
def build(Mtot,mtab,nT=500):
    X0,MT=mtab;Mx=lambda x:Mtot*np.interp(x,X0,MT)
    E0=0.5*u*u-Mtot;t0=age_out(1.0,Mtot,E0)
    XL=np.r_[np.geomspace(2e-3,1.5,500),np.geomspace(1.5,40,200)[1:]]
    ok=np.array([Mx(x)==0 or age_out(x,Mx(x),-Mx(x)/x*(1-1e-12))>=t0 for x in XL])   # shells still on the outgoing branch today
    nbad=int((~ok).sum())
    if nbad and ok[:np.argmax(ok)].any()==False: XL=XL[np.argmax(ok):]
    elif nbad: raise RuntimeError('non-contiguous turnaround shells')
    E=np.empty_like(XL);Mv=np.array([Mx(x) for x in XL])
    for i,(x,m) in enumerate(zip(XL,Mv)):
        if m==0: E[i]=0.5*(x/t0)**2;continue
        lo=-m/x*(1-1e-12);
        E[i]=brentq(lambda e:age_out(x,m,e)-t0,lo,50*(x/t0)**2+10*m/x+10,xtol=1e-15,rtol=1e-13)
    T=np.geomspace(2e-3,t0,nT)
    W=np.sign(E)*np.sqrt(2*np.abs(E))
    G=np.array([[V.Rfun(t,w,m) for w,m in zip(W,Mv)] for t in T])
    cross=bool(np.any(np.diff(G,axis=1)<=0))
    S=RectBivariateSpline(T,XL,G,kx=5,ky=5)
    tc=brentq(lambda t:S(t,1.0)[0,0]-R,T[0]*1.01,t0)
    return dict(n_dropped_inner=nbad,xmin=float(XL[0]),t0=t0,tc=tc,T=T,X=XL,S=S,E=E,M=Mv,W=W,shell_crossing=cross,Mtot=Mtot,G=G)
def trace(F,mode,psis):
    S,T,W,t0=F['S'],F['T'],F['X'],F['t0'];wo=1.0
    def AB(t,w):
        Rr=S(t,w)[0,0];Rp=S(t,w,dy=1)[0,0];Rt=S(t,w,dx=1)[0,0];Rpp=S(t,w,dy=2)[0,0];Rtp=S(t,w,dx=1,dy=1)[0,0]
        return Rr,Rp*Rp,2*Rp*Rtp,2*Rp*Rpp,Rr*Rr,2*Rr*Rt,2*Rr*Rp
    def rhs(l,y):
        t,w,ph,pt,pw,pp=y;Rr,A,At,Aw,B,Bt,Bw=AB(t,w)
        return [-pt,pw/A,pp/B,0.5*At/A**2*pw**2+0.5*Bt/B**2*pp**2,0.5*Aw/A**2*pw**2+0.5*Bw/B**2*pp**2,0]
    if mode=='own': hit=lambda l,y:AB(y[0],y[1])[0]-R
    else: hit=lambda l,y:y[0]-F['tc']
    hit.terminal=True
    edge=lambda l,y:min(y[0]-T[0]*1.01,y[1]-W[0]-1e-3,W[-1]-1e-2-y[1]);edge.terminal=True
    A0=AB(t0,wo)[1];B0=AB(t0,wo)[4];out=[]
    for ps in psis:
        pr=np.radians(ps);y0=[t0,wo,0.,-1.,-np.sqrt(A0)*np.cos(pr),-np.sqrt(B0)*np.sin(pr)]
        s=solve_ivp(rhs,(0,-200),y0,events=[hit,edge],rtol=1e-9,atol=1e-11)
        if s.t_events[0].size:
            y=s.y_events[0][0];out.append((ps,-y[3],y[0],t0-y[0],y[1],y[2],S(y[0],y[1])[0,0]))
        else: out.append((ps,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan))
    return np.array(out)
if __name__=='__main__':
    tag,kind,frac,mode,step,outf=sys.argv[1:7];frac=float(frac);step=float(step)
    LD=18159.848409348455/39.10676196793444;alpha=0.30000009809958916;DR=42/14.30048287260589;Rb=1.5*(1+1/DR)
    mt=mfrac_table(LD,alpha,Rb,kind);F=build(frac*0.5*u*u,mt)
    Gy=42*3.26156
    print(tag,'t0 Gyr',F['t0']*Gy,'our crossing Gyr',F['tc']*Gy,'shell crossing on grid',F['shell_crossing'],flush=True)
    o=trace(F,mode,np.arange(0,180.001,step))
    np.savez(outf,o=o,t0=F['t0'],tc=F['tc'],Mtot=F['Mtot'],shell_crossing=F['shell_crossing'],LD=LD,alpha=alpha,Rb=Rb,kind=kind,frac=frac,mode=mode)
    print(tag,'resolved',int(np.isfinite(o[:,1]).sum()),'/',len(o),'1+z',np.nanmin(o[:,1]),np.nanmax(o[:,1]),flush=True)

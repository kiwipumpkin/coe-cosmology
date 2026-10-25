"""Axially symmetric quasispherical Szekeres (S(r) only, P=Q=0; Buckley & Schlegel 2019 eqs 1,2,6,7), dipole axis = CoE axis.
Plane containing the axis (polar angle psi from the axis, psi=0 toward the observer / away from CoE):
 ds^2 = -dt^2 + (Phi' + Phi s cos psi)^2/(1+2E) dr^2 + Phi^2 (dpsi - s sin psi dr)^2,  s(r) = S'/S.
Built on the LTB toy (valley/toy.py) for Phi(r,t), M(r), E(r)."""
import sys,json,pickle,warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/valley'); sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/ltb'); sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr')
import numpy as np, toy
from ltbcore import _ev
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, fsolve
from grcore import relsum
def smooth_bump(r,a,b,w):
    x=np.clip((r-a)/w,0,1); y=np.clip((b-r)/w,0,1); f=lambda z: z*z*z*(10-15*z+6*z*z); return f(x)*f(y)
class SZ:
    def __init__(s,L,C,a,b,w=0.15,sig0=None):
        s.L=L; rg=L.rg; M=np.array([L.Mf(r) for r in rg]); Mp=np.gradient(M,rg)
        s.bump=smooth_bump(rg,a,b,w); inside=s.bump>0.999
        tt=L.tg[(L.tg>=L.t1-0.5)&(L.tg<=L.t1+4)]; Pr=L.Rs(rg,tt,dx=1); P=L.Rs(rg,tt)
        s.limM=Mp/(3*M); s.limP=(Pr/P).min(1)       # no-shell-crossing limits on |S'/S|
        s.sig0=sig0 if sig0 is not None else float(min(s.limM[inside].min(),s.limP[inside].min()))
        sg=C*s.sig0*s.bump; s.sg=CubicSpline(rg,sg); s.C=C
        s.ratio=float(np.max(np.abs(sg)/np.minimum(np.maximum(s.limM,1e-30),s.limP)))   # must be < 1
        s.cint=CubicSpline(rg,cumulative_trapezoid(rg*sg,rg,initial=0.0))           # shell-centre offset at t1 (Phi=r there)
        s.t1=L.t1; s.rg=rg
    def F(s,t,r):
        L=s.L; P=L.Rs.ev(r,t); Pr=_ev(L.Rr,r,t); Pt=_ev(L.Rt_,r,t); Prt=_ev(L.Rrt,r,t); Prr=_ev(L.Rrr,r,t)
        E=L.Es(r); Ep=L.Es(r,1); return P,Pr,Pt,Prt,Prr,1+2*E,2*Ep,s.sg(r),s.sg(r,1)
    def metric(s,t,r,ps,d=False):
        P,Pr,Pt,Prt,Prr,K,Kr,g,gp=s.F(t,r); c,sn=np.cos(ps),np.sin(ps); A=Pr+P*g*c
        G=np.array([[A*A/K+P*P*g*g*sn*sn, -P*P*g*sn],[-P*P*g*sn, P*P]])
        if not d: return G,A,K,P
        At=Prt+Pt*g*c; Ar=Prr+Pr*g*c+P*gp*c; Ap=-P*g*sn
        Gt=np.array([[2*A*At/K+2*P*Pt*g*g*sn*sn, -2*P*Pt*g*sn],[-2*P*Pt*g*sn, 2*P*Pt]])
        Gr=np.array([[2*A*Ar/K-A*A*Kr/K**2+2*P*Pr*g*g*sn*sn+2*P*P*g*gp*sn*sn, -(2*P*Pr*g+P*P*gp)*sn],[-(2*P*Pr*g+P*P*gp)*sn, 2*P*Pr]])
        Gp=np.array([[2*A*Ap/K+2*P*P*g*g*sn*c, -P*P*g*c],[-P*P*g*c, 0.0]])
        return G,Gt,Gr,Gp
    def rhs(s,t,y):
        r,ps,vr,vp,_=y; G,Gt,Gr,Gp=s.metric(t,r,ps,True); v=np.array([vr,vp]); Gi=np.linalg.inv(G); dG=[Gr,Gp]
        Gtt=0.5*v@Gt@v
        Gam=np.zeros((2,2,2))   # Gamma^i_jk spatial
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    Gam[i,j,k]=0.5*sum(Gi[i,l]*(dG[j][l,k]+dG[k][l,j]-dG[l][j,k]) for l in range(2))
        a=-(Gi@Gt@v)-np.einsum('ijk,j,k->i',Gam,v,v)+v*Gtt
        return [vr,vp,a[0],a[1],-Gtt]
    def frame(s,t,r,ps,wr,wp):   # local radial/tangential (w.r.t. CoE) velocity measured by the comoving observer -> coordinate velocity
        G,A,K,P=s.metric(t,r,ps); N=np.sqrt(K)/A; g=s.sg(r); return np.array([wr*N, wr*N*g*np.sin(ps)+wp/P])
    def kdot(s,t,r,ps,va,lna,vb,lnb,null_b):
        G,_,_,_=s.metric(t,r,ps); q=-1+va@G@vb
        if null_b: return np.exp(lna+lnb)*q
        return np.exp(lna)*q/np.sqrt(1-vb@G@vb)
    # Euclidean "shifted shells" picture at t1 (exact LTB labels when the dipole is zero)
    def X(s,r,ps): return np.array([r*np.cos(ps)+s.cint(r), r*np.sin(ps)])
    def inv(s,Xv):
        r=np.hypot(*Xv)
        for _ in range(60): x0=Xv[0]-s.cint(r); r=np.hypot(x0,Xv[1])
        return r,np.arctan2(Xv[1],x0)
def edge(S):
    ev=lambda t,y: min(y[0]-S.rg[0]-0.2, S.rg[-1]-1-y[0]); ev.terminal=True; return ev
def wl(S,X0,w,ts,T):
    r,ps=S.inv(X0); wr=w[0]*np.cos(ps)+w[1]*np.sin(ps); wp=-w[0]*np.sin(ps)+w[1]*np.cos(ps); v=S.frame(ts,r,ps,wr,wp)
    so=solve_ivp(S.rhs,(ts,T),[r,ps,*v,0.0],method="DOP853",dense_output=True,rtol=1e-10,atol=1e-12,events=edge(S)); te=so.t[-1]
    return lambda t: None if (t<ts-1e-12 or t>te+1e-12) else so.sol(t)
def photon(S,ro,po,b,T,tl):
    v=S.frame(T,ro,po,np.cos(b),np.sin(b)); return solve_ivp(S.rhs,(T,tl),[ro,po,*v,0.0],method="DOP853",dense_output=True,rtol=1e-10,atol=1e-12,events=edge(S)),v
def ring(S,v0,T,D):
    g=1/np.sqrt(1-v0*v0); X0=S.X(D,0.0); ems=[]
    for t in np.radians(np.arange(0,361,1.0)):
        e=np.array([np.cos(t),np.sin(t)]); ts=S.t1+g*v0*(1+e[0]); w=relsum(v0,1.0,e); w/=np.linalg.norm(w)
        ems.append((wl(S,X0+np.array([g*(e[0]+v0*v0),e[1]]),w,ts,T),ts,g*(1+v0*e[0])))
    ob=wl(S,X0+np.array([g*v0*v0,0]),np.array([v0,0.0]),S.t1+g*v0,T)
    return ems,ob(T)
def observe(S,y,ems,T):
    ro,po,vro,vpo,_=y; out=[]; prev=None; Xo=S.X(ro,po)
    for k,(em,ts,l0) in enumerate(ems):
        hi=T-1e-6
        while em(hi) is None: hi=ts+(hi-ts)*0.99
        e0=em(0.5*(ts+hi)); dX=S.X(e0[0],e0[1])-Xo; b0=np.arctan2(-(dX@[-np.sin(po),np.cos(po)]),-(dX@[np.cos(po),np.sin(po)]))
        g=lambda s_: np.linalg.norm(S.X(*em(s_)[:2])-Xo)-(T-s_)
        try: s0=brentq(g,ts,hi)
        except Exception: s0=ts
        def res(z):
            b,s_=z
            if not (ts<=s_<=T): return [10+abs(s_),10]
            so,_=photon(S,ro,po,b,T,ts-0.2)
            if so.t[-1]>s_: return [5,5]
            e=em(s_)
            if e is None: return [7,7]
            yy=so.sol(s_); return [yy[0]-e[0], e[0]*(yy[1]-e[1])]
        sol=None
        for gs in ([b0,s0],prev):
            if gs is None: continue
            z=fsolve(res,gs,xtol=1e-12)
            if np.max(np.abs(res(z)))<1e-8: sol=z; break
        if sol is None: out.append(np.nan); continue
        prev=sol; b,s_=sol; so,vph=photon(S,ro,po,b,T,ts-0.2); yy=so.sol(s_); e=em(s_)
        num=S.kdot(T,ro,po,vph,0.0,np.array([vro,vpo]),0.0,False)
        den=S.kdot(s_,yy[0],yy[1],yy[2:4],yy[4],e[2:4],e[4]+np.log(l0),True)
        out.append(num/den)
    return np.array(out)
def migkas_af(S,D,TH):   # expansion-rate pattern seen from Earth over the vision sphere (radius 1 = RV) at t1, shifted-shell picture
    L=S.L; Xe=S.X(D,0.0); cd=CubicSpline(S.rg,cumulative_trapezoid(np.array([L.Rs.ev(r,S.t1,dy=1) for r in S.rg])*S.sg(S.rg),S.rg,initial=0.0))
    V=lambda r,ps: L.Rs.ev(r,S.t1,dy=1)*np.array([np.cos(ps),np.sin(ps)])+np.array([cd(r),0.0])
    Ve=V(D,0.0); af=[]
    for th in np.radians(TH):
        r,ps=S.inv(Xe+np.array([np.cos(th),np.sin(th)])); af.append(np.linalg.norm(V(r,ps)-Ve))
    af=np.array(af); return af/af.mean()

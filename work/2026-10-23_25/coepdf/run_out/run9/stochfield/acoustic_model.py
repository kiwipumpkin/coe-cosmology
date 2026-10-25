# Acoustic-feature stochastic field: exact projection for l<=LX, Limber for l>LX (moving-lens term dominant; SW/ISW scalar kept only in exact part).
import json,numpy as np
from fwd import Fwd,T0,Gc,rho
LX=50;LMAX=2500
class AFwd(Fwd):
    def Pd(s,p,k=None):
        k=s.k if k is None else k
        return 10**p['lgA']*k**p['n']*(1+p['Aw']*np.cos(k*p['rs']+p['phi']))*np.exp(-(k/p['kD'])**2)
EX=AFwd(lmax=LX,kmax=200,nk=500)
lL=np.arange(LX+1,LMAX+1)
def limber(p,l=lL,n=600):
    a=max(0.01,p['chiw']-6*p['sw']);b=min(EX.R,p['chiw']+6*p['sw'])
    if b<=a+1e-6:return np.zeros(len(l))
    ch=np.linspace(a,b,n);KM=np.interp(ch,EX.ch,EX.KM)
    W=np.exp(-0.5*((ch-p['chiw'])/p['sw'])**2);Wd=-W*(ch-p['chiw'])/p['sw']**2
    k=(l[:,None]+0.5)/ch[None];Pp=(4*np.pi*Gc*rho)**2*EX.Pd(p,k)/k**4
    I=(KM*W)**2*k**2/3*Pp/ch**2   # Limber: only transverse k; <(k.u)^2> over sky = k^2/3; W' term is part of the radial (suppressed) modes
    return np.trapezoid(I,ch,axis=1)*l*(l+1)/(2*np.pi)*T0**2
def Dfield(p):
    t=EX.terms(p);return np.r_[0,0,t['scalar']+t['mlens'],limber(p)]   # index = l, 0..2500
Dfix=np.zeros(LMAX+1);Dfix[2:31]=json.load(open('/workspace/redweb2/stochfield/fit_l30_bg.json'))['Dfix']
def Dmodel(p):Df=Dfield(p);return Df+Dfix,Df
from scipy.special import spherical_jn
def exact_nodes(p,nodes=(51,60,75,90,110,130,150,175,200,240,280,330,400)):
    """large-l-limit exact projection (validated vs Fwd exact to <2% at l 20-50) at nodes; returns ratio exact/limber"""
    ch,KM,dch=EX.ch,EX.KM,EX.dch;W=np.exp(-0.5*((ch-p['chiw'])/p['sw'])**2);Wd=-W*(ch-p['chiw'])/p['sw']**2;m=W>1e-4;ch,KM,W,Wd=ch[m],KM[m],W[m],Wd[m]
    r=[]
    for L in nodes:
        k=np.linspace(L/ch.max()*0.7,min((L+60)/ch.min()*2.5,4000),2500);x=k[:,None]*ch[None]
        A=spherical_jn(L,x)@(KM*W)*dch;Bp=spherical_jn(L+1,x)@(KM*Wd)*dch;Bm=spherical_jn(L-1,x)@(KM*Wd)*dch
        acc=k**2*A**2/3+(Bp**2+Bm**2)/6+k*A*(Bp-Bm)/3;Pp=(4*np.pi*Gc*rho)**2*EX.Pd(p,k)/k**4
        C=2/np.pi*np.trapezoid(k**2*Pp*acc,k)*L*(L+1)/(2*np.pi)*T0**2;r.append(C/limber(p,np.array([L]))[0])
    return np.array(nodes),np.array(r)
def corr_curve(p):
    n,r=exact_nodes(p);c=np.ones(LMAX+1);ll=np.arange(LMAX+1);sel=(ll>LX)&(ll<=400);c[sel]=np.exp(np.interp(ll[sel],n,np.log(r)));return c,n,r

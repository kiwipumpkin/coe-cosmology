# Semi-analytic forward model: ensemble C_l of dT/T from a Gaussian random potential field (envelope W(chi)) riding the CoE flow.
# Kernels copied from los_flow_R/flowR.py (grow variant, full CoE radial velocity in moving lens), evaluated at x=1 (CoE-radius of observer).
import numpy as np,healpy as hp
from scipy.special import spherical_jn
from scipy.interpolate import RegularGridInterpolator as RGI
c=2.998e5;Gc=4.30e-12;rho=42152249999.99999*1e9;T0=2.7255e6
T=np.load('/workspace/redweb2/realweb/los_flow_R/tab.npz');L1=np.load('/workspace/redweb2/realweb/los_flow/ltb.npz')
lTs=np.log(T['Ts']);Stf=RGI((lTs,T['xs']),T['St']);ix1=np.argmin(abs(T['xs']-1));a_=T['S'][:,ix1];eta=T['eta']
A=lambda t:np.interp(np.log(t),np.log(L1['T']),L1['Agrow'])
def t_of_chi(ch,R):return np.exp(np.interp(-ch/R*eta[0],-eta,lTs))
def kernels(R,dch=0.004):
    ch=np.arange(0,R,dch)+dch/2;te=t_of_chi(np.r_[ch-dch/2,R],R);tt=t_of_chi(ch,R);dt=-np.diff(te)
    ai=np.interp(np.log(tt),lTs,a_);Am=0.5*(A(te[:-1])+A(te[1:]))
    KI=(2/c**2)*(A(te[:-1])-A(te[1:]))/dch                           # per (km/s)^2 per Gpc
    Rd=Stf(np.c_[np.log(tt),np.ones_like(tt)])
    KM=(2/c**2)*(-(Am/ai*dt)*42*Rd)/dch                               # multiplies u.grad Phi0 [(km/s)^2/Gpc]
    SW=A(t_of_chi(R-1e-9,R))/(3*c**2)
    return ch,KI,KM,SW
class Fwd:
    def __init__(s,R=4.0,lmax=100,nk=360,kmin=1e-2,kmax=120,dch=0.004):
        s.R=R;s.l=np.arange(2,lmax+1);s.lx=np.arange(0,lmax+2);s.k=np.geomspace(kmin,kmax,nk);s.lnk=np.log(s.k)
        s.ch,s.KI,s.KM,s.SW=kernels(R,dch);s.dch=dch
        s.J=np.empty((len(s.lx),nk,len(s.ch)),np.float32);s.JR=np.empty((len(s.lx),nk))
        for i,L in enumerate(s.lx):s.J[i]=spherical_jn(L,s.k[:,None]*s.ch[None]);s.JR[i]=spherical_jn(L,s.k*R)
    def W(s,p,der=False):  # radial envelope on potential: Gaussian in distance, centre chiw, width sw (Gpc)
        w=np.exp(-0.5*((s.ch-p['chiw'])/p['sw'])**2)
        if der:return -w*(s.ch-p['chiw'])/p['sw']**2
        return w,np.exp(-0.5*((s.R-p['chiw'])/p['sw'])**2)
    def Pd(s,p,k=None):  # density power spectrum P_delta(k) [Gpc^3]: A (k/1)^n exp(-(k/kc)^2)
        k=s.k if k is None else k;return 10**p['lgA']*k**p['n']*np.exp(-(k/p['kc'])**2)
    def Pphi(s,p):return (4*np.pi*Gc*rho)**2*s.Pd(p)/s.k**4
    def terms(s,p):
        w,wR=s.W(p);wd=s.W(p,True);Is=(np.einsum('lkc,c->lk',s.J,(s.KI*w).astype(np.float32))*s.dch+s.SW*wR*s.JR)[2:-1]
        Am=np.einsum('lkc,c->lk',s.J,(s.KM*w).astype(np.float32))*s.dch;Bm=np.einsum('lkc,c->lk',s.J,(s.KM*wd).astype(np.float32))*s.dch
        P=s.Pphi(p);f=2/np.pi*s.k**3*P   # integrand per dlnk
        Cs=np.trapezoid(f*Is**2,s.lnk,axis=1);Cm=np.zeros(len(s.l));k=s.k
        # moving lens  int K_M u.grad(W phi) = W u.grad phi + phi W'(n.u): exact m-averaged ensemble C_l via cos(theta) Gaunt couplings (u = z)
        if not hasattr(s,'cf'):
            a=lambda l,m:np.sqrt(((l+1)**2-m**2)/((2*l+1)*(2*l+3)));b=lambda l,m:np.sqrt(np.maximum(l**2-m**2,0)/((2*l-1)*(2*l+1)))
            cf=[]
            for l in s.l:
                m=np.arange(-l,l+1);mi=np.abs(m)<=l-1
                cf.append([np.sum(a(l,m)**2+b(l,m)**2),np.sum(b(l+1,m)**2),np.sum(a(l-1,m[mi])**2),np.sum(a(l,m)*b(l+1,m)),np.sum(b(l,m[mi])*a(l-1,m[mi]))])
            s.cf=np.array(cf)/(2*s.l[:,None]+1)
        L=s.l;A0=Am[L];Bp=Bm[L+1];Bn=Bm[L-1];c=s.cf
        acc=k**2*A0**2*c[:,:1]+Bp**2*c[:,1:2]+Bn**2*c[:,2:3]+2*k*A0*Bp*c[:,3:4]-2*k*A0*Bn*c[:,4:5]
        Cm=np.trapezoid(f*acc,s.lnk,axis=1)
        dl=s.l*(s.l+1)/(2*np.pi)*T0**2;return dict(scalar=Cs*dl,mlens=Cm*dl)
    def Dl(s,p):t=s.terms(p);return t['scalar']+t['mlens']
    def sig_delta(s,p,Rs=0.078):  # rms density contrast smoothed with Gaussian of Rs Gpc
        kk=np.geomspace(1e-3,1e3,4000);return float(np.sqrt(np.trapezoid(kk**3*s.Pd(p,kk)*np.exp(-(kk*Rs)**2)/(2*np.pi**2),np.log(kk))))

# Exact line-of-sight projection (large-l form, validated vs Fwd exact at l=30-50 to <0.2%) with k grid resolving j_l oscillations
import numpy as np
from scipy.special import spherical_jn
import acoustic_model as M
def node(args):
    p,L,dk=args;EX=M.EX
    ch,KM,dch=EX.ch,EX.KM,EX.dch;W=np.exp(-0.5*((ch-p['chiw'])/p['sw'])**2);Wd=-W*(ch-p['chiw'])/p['sw']**2;m=W>1e-4;ch,KM,W,Wd=ch[m],KM[m],W[m],Wd[m]
    k0=L/ch.max()*0.7;k1=min((L+60)/ch.min()*2.5,4000);k=np.linspace(k0,k1,max(2500,int((k1-k0)/dk)))
    A=np.zeros(len(k));Bp=np.zeros(len(k));Bm=np.zeros(len(k))
    for i in range(0,len(k),2000):
        x=k[i:i+2000,None]*ch[None];A[i:i+2000]=spherical_jn(L,x)@(KM*W)*dch;Bp[i:i+2000]=spherical_jn(L+1,x)@(KM*Wd)*dch;Bm[i:i+2000]=spherical_jn(L-1,x)@(KM*Wd)*dch
    acc=k**2*A**2/3+(Bp**2+Bm**2)/6+k*A*(Bp-Bm)/3;Pp=(4*np.pi*M.Gc*M.rho)**2*M.EX.Pd(p,k)/k**4
    C=2/np.pi*np.trapezoid(k**2*Pp*acc,k)*L*(L+1)/(2*np.pi)*M.T0**2;return C/M.limber(p,np.array([L]))[0]

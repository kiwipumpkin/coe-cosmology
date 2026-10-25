import numpy as np,json
from acoustic_model import *
from scipy.special import spherical_jn
out=[]
PS=[dict(lgA=-9,n=0.5,Aw=0.5,rs=0.084,phi=0,kD=300,chiw=c,sw=w) for c,w in ((3.95,0.03),(2.5,0.5),(4.0,0.1),(3.7,0.1),(4.0,0.2))]
import sys;PS=PS[2:] if len(sys.argv)>1 else PS
for p in PS:
    D=Dmodel(p)[1];ch,KM,dch=EX.ch,EX.KM,EX.dch;W=np.exp(-0.5*((ch-p['chiw'])/p['sw'])**2);Wd=-W*(ch-p['chiw'])/p['sw']**2;m=W>1e-4;ch,KM,W,Wd=ch[m],KM[m],W[m],Wd[m]
    for L in ((51,70,100,150,220,300,540) if len(sys.argv)>1 else (20,40,50,70,100,150,220,300,540,800,1500)):
        k=np.linspace(L/ch.max()*0.7,min((L+60)/ch.min()*2.5,4000),6000);x=k[:,None]*ch[None]
        A=spherical_jn(L,x)@(KM*W)*dch;Bp=spherical_jn(L+1,x)@(KM*Wd)*dch;Bm=spherical_jn(L-1,x)@(KM*Wd)*dch
        acc=k**2*A**2/3+(Bp**2+Bm**2)/6+k*A*(Bp-Bm)/3      # large-l limit of the exact Gaunt-coupled m-average
        Pp=(4*np.pi*Gc*rho)**2*EX.Pd(p,k)/k**4;C=2/np.pi*np.trapezoid(k**2*Pp*acc,k)*L*(L+1)/(2*np.pi)*T0**2
        Ll=limber(p,np.array([L]))[0];out.append(dict(chiw=p['chiw'],sw=p['sw'],l=L,exact=C,limber=Ll));print(p['sw'],L,'exact %.4g limber %.4g ratio %.3f'%(C,Ll,Ll/C),('fwd-exact %.4g'%D[L] if L<=50 else ''),flush=True)
json.dump(out,open('acoustic_limber_check'+('2' if len(sys.argv)>1 else '')+'.json','w'),indent=1)

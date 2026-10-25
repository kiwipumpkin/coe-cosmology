import sys,numpy as np,base as V
u=float(sys.argv[1]);f=float(sys.argv[2]);V.uD=u
import v11
from scipy.optimize import brentq
M=f*0.5*u**2;R=0.5;F=v11.field(M,R);wo,t0,T,W,S=F;o,_=v11.trace(F,R,np.arange(0,181,1.0))
re=np.array([V.Rfun(a,w,M) for a,w in zip(o[:,2],o[:,4])]);rt=np.array([V.Rfun(t0,w,M) for w in o[:,4]])
xo=V.Rfun(t0,wo,M);tc=brentq(lambda t:V.Rfun(t,wo,M)-R,T[0]*1.01,t0);xc=V.Rfun(tc,wo,M)
np.savez(sys.argv[3],o=o,re=re,rt=rt,phi=o[:,5],xo=xo,xc=xc,t0=t0,tc=tc);print(sys.argv[3],np.isfinite(o[:,1]).sum())

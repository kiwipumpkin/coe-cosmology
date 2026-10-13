import numpy as np,json,sys
from common import *
from fit_on_web import spec
K=seq_core.unit_pattern_fast; tag=sys.argv[1]
f=np.load(f"final_{tag}.npz"); X0=f["Xp"]; r0=f["rp"]; isV=f["kp"]=="void"; M0=f["Mp"]; n=len(X0)
C=M0[:,None]*K(X0)/FLOOR; b=(R-f["w"])/FLOOR
U,S,Vt=np.linalg.svd(C,full_matrices=False); print("sing values (max, 100th, 300th, min):",S[0],S[99],S[299],S[-1])
x,*_=np.linalg.lstsq(C.T,b,rcond=None); res=b-C.T@x; print("unbounded LS: left %.3g, max|a-1| %.3g, median|da| %.3g"%(rms(res),np.abs(x).max(),np.median(np.abs(x))))
# component of b outside the span of the 360-dim kernel family? constant offset check
print("b mean",b.mean(),"col means", np.abs(C.mean(1)).max()/np.abs(C).max())
dr=np.sqrt((np.hypot(X0[:,0],X0[:,1])-RMIN)**2+X0[:,2]**2); near=(dr-r0)<3
x2,*_=np.linalg.lstsq(C[near].T,b,rcond=None); print("near-ring only unbounded: left %.3g max|da| %.3g"%(rms(b-C[near].T@x2),np.abs(x2).max()))

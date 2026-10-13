import numpy as np,json,sys
from scipy.optimize import lsq_linear
from common import *
from fit_on_web import spec
K=seq_core.unit_pattern_fast
tag=sys.argv[1]
f=np.load(f"final_{tag}.npz"); X0=f["Xp"]; r0=f["rp"]; isV=f["kp"]=="void"; M0=f["Mp"]; n=len(X0)
d=json.load(open(f"minedit2/{tag}.json")); i=np.array(d["idx"]); XS=X0.copy(); AS=np.ones(n); XS[i]=d["X_new"]; AS[i]=d["a_new"]
rho0=seq_core.rho_gpc(X0); m=np.where(isV,M0*seq_core.rho_gpc(XS)/rho0,M0)
C=m[:,None]*K(XS); W=(AS[:,None]*C).sum(0); b=(R-W)/FLOOR
dr=np.sqrt((np.hypot(XS[:,0],XS[:,1])-RMIN)**2+XS[:,2]**2)
for nm,amax in (("shrink/delete only",AS),("grow to x2",np.maximum(AS,2)),("grow to x10",np.maximum(AS,10))):
    lo=-AS; hi=amax-AS; act=hi-lo>1e-9
    from scipy.sparse.linalg import aslinearoperator
    s=lsq_linear(aslinearoperator(np.ascontiguousarray(C[act].T/FLOOR)),b,bounds=(lo[act],hi[act]),method="trf",lsmr_tol="auto",max_iter=300,verbose=0)
    x=np.zeros(n); x[act]=s.x; res=(b-C.T@x/FLOOR)*FLOOR; sp=spec(res)
    ch=np.abs(x)>1e-3
    print(tag,nm,"left %.2f"%(rms(res)/FLOOR),"m13 %.2f"%sp["rms_13up"],"changed",ch.sum(),"near-ring changed",(ch&((dr-r0)<3)).sum())

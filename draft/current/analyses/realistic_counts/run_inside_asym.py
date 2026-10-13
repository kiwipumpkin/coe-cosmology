"""Inside-sphere point fits with the no-mirror rule of ../solution_space/asym (box rule, 20 deg / 15% / 25%, 10 deg off plane)."""
import os, sys, json, numpy as np
os.environ["OMP_NUM_THREADS"]="1"
CS=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,os.path.join(CS,"inside_visible")); sys.path.insert(0,os.path.join(CS,"solution_space","asym"))
from inside_core import *
import asym_core as AC
from scipy.optimize import least_squares
from multiprocessing import Pool
SC=1/(FLOOR*np.sqrt(360))
def job(a):
    kinds,seed=a; m=InModel(kinds); rng=np.random.default_rng(seed); n=m.n
    p=np.zeros((n,4)); p[:,0]=rng.uniform(0,1,n); p[:,1]=rng.uniform(-np.pi,np.pi,n); p[:,2]=np.arcsin(rng.uniform(-1,1,n)); p[:,3]=rng.uniform(0,1,n); p=p.ravel()
    lam=30.0
    for it in range(5):
        def f(q):
            X,M,Rv,Rc=m.decode(q); return np.concatenate([(R-m.pattern(q))*SC,lam*m.violations(q),lam*AC.asym_viol(m.k,X,M,Rv)])
        p=least_squares(f,p,bounds=(m.lo,m.hi),method="trf",x_scale="jac",max_nfev=60*len(p),xtol=1e-8,ftol=1e-8).x
        X,M,Rv,Rc=m.decode(p)
        if m.violations(p).max(initial=0)<=1e-3 and AC.asym_viol(m.k,X,M,Rv).max(initial=0)<=1e-3: break
        lam*=10
    X,M,Rv,Rc=m.decode(p); nm,npl,_,_=AC.hard_check(m.k,X,M,Rv)
    return dict(kinds=kinds,seed=seed,left=float(m.leftover(p)),viol=float(m.violations(p).max(initial=0)),mirror_matches=nm,on_plane=npl,X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),Rc=Rc.tolist())
if __name__=="__main__":
    ks=["CCCCCCC","CCCCCCCC","CCCCCCVV","CCCCVVVV","CCCCCCCCCC"]
    with Pool(4) as pool: out=pool.map(job,[(k,100+i) for k in ks for i in range(10)])
    json.dump(out,open("inside_asym.json","w"))
    for k in ks:
        ok=[r for r in out if r["kinds"]==k and r["viol"]<=1e-3 and r["mirror_matches"]==0 and r["on_plane"]==0]
        print(k,"valid",len(ok),"best",round(min([r["left"] for r in ok],default=np.nan),3))

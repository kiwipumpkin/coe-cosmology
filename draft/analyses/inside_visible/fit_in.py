"""Multi-start fits, objects inside the visible sphere. usage: python3 fit_in.py TAG NSTARTS MHI KINDS1 KINDS2 ... -> raw_TAG/KINDS.json"""
import os, sys, json, time, numpy as np
os.environ["OMP_NUM_THREADS"]="1"
from multiprocessing import Pool
from scipy.optimize import least_squares
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from inside_core import *
SC=1/(FLOOR*np.sqrt(360))
def fit_one(m,p0):
    lam=30.0; p=p0
    for it in range(4):
        f=lambda q: np.concatenate([(R-m.pattern(q))*SC, lam*m.violations(q)])
        r=least_squares(f,p,bounds=(m.lo,m.hi),method="trf",x_scale="jac",max_nfev=60*len(p),xtol=1e-8,ftol=1e-8)
        p=r.x
        if m.violations(p).max(initial=0)<=1e-3: break
        lam*=10
    return p
def job(a):
    kinds,mhi,seed=a; m=InModel(kinds,mhi,os.environ.get("COMP","local")); rng=np.random.default_rng(seed); n=m.n
    p0=np.zeros((n,4)); p0[:,0]=rng.uniform(0,1,n); p0[:,1]=rng.uniform(-np.pi,np.pi,n); p0[:,2]=np.arcsin(rng.uniform(-1,1,n)); p0[:,3]=rng.uniform(0,1,n)
    t=time.time(); p=fit_one(m,p0.ravel()); X,M,Rv,Rc=m.decode(p); w=m.pattern(p)
    return dict(kinds=kinds,seed=seed,comp=m.comp,left=float(rms(R-w)/FLOOR),left_localcomp=float(m.leftover(p,"local")),left_surroundcomp=float(m.leftover(p,"surround")),left_nocomp=float(m.leftover(p,"none")),
        viol=float(m.violations(p).max(initial=0)),p=p.tolist(),X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),Rc=Rc.tolist(),amps=(amps(R-w)/FLOOR).tolist(),secs=time.time()-t)
if __name__=="__main__":
    tag,ns,mhi=sys.argv[1],int(sys.argv[2]),float(sys.argv[3]); ks=sys.argv[4:]; os.makedirs(f"raw_{tag}",exist_ok=True)
    tasks=[(k,mhi,1000*len(k)+i) for k in ks for i in range(ns) if not os.path.exists(f"raw_{tag}/{k}.json")]
    out={}
    with Pool(int(os.environ.get("NP","6"))) as pool:
        for r in pool.imap_unordered(job,tasks):
            out.setdefault(r["kinds"],[]).append(r)
            if len(out[r["kinds"]])==ns:
                runs=sorted(out[r["kinds"]],key=lambda q:q["seed"]); json.dump(dict(kinds=r["kinds"],mhi=mhi,nstarts=ns,runs=runs),open(f"raw_{tag}/{r['kinds']}.json","w"))
                L=sorted(q["left"] for q in runs if q["viol"]<=1e-3)
                print(tag,r["kinds"],f"best={L[0] if L else float('nan'):.2f} med={np.median(L) if L else float('nan'):.1f} feas={len(L)}/{ns} t/fit={np.mean([q['secs'] for q in runs]):.1f}s",flush=True)

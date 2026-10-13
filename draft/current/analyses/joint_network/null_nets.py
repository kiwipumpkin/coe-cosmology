"""Random 'LCDM-like' networks under the same rules, to compare with a fitted layout.
For a fitted layout (level, mirror, N): K networks with the SAME bodies (masses, void radii, depths shuffled within kind),
placed isotropically and homogeneously (uniform in volume) between the fit's nearest and farthest body distance,
then relaxed with the rule penalties only (no CMB term), using the fit's own rule level and weights.
usage: python3 null_nets.py NPROC K level sym N  -> null/<level>_<sym>_N<N>.json"""
import os,sys,json,time,numpy as np
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"): os.environ[v]="1"
from multiprocessing import Pool
from scipy.optimize import least_squares
from jn_fast import JF
from jn_core import *
def best(level,sym,N):
    r=[x for x in json.load(open(f"raw/{level}_{sym}_N{N}.json")) if x["ok"]]
    lam=LEVELS[level]["lam"]; return min(r,key=lambda x:x["left"]**2+(lam*x["soft_rms"])**2*len(x["p"]))
def relax(m,p,lam,lamh=30.0):
    for it in range(3):
        w=m.wvec(lam,lamh)
        r=least_squares(lambda q:m.fun(q,w)[360:],p,jac=lambda q:m.jac(q,w)[360:],bounds=(m.lo,m.hi),method="trf",x_scale="jac",max_nfev=300,xtol=1e-9,ftol=1e-9)
        p=r.x; h,a,s=m.split(m.decode(p))
        if h.max(initial=0)<=1e-3 and a.max(initial=0)<=1e-3: break
        lamh*=10
    return p
def one(a):
    level,mir,kinds,X,M,Rv,q,dmin,dmax,seed=a; rng=np.random.default_rng(seed); m=JF(kinds,level,mir)
    k=np.array(list(kinds)); n=len(k); X=np.array(X); M=np.abs(np.array(M)); Rv=np.array(Rv); q=np.array(q)
    for kk in "CV":
        i=np.where(k==kk)[0]; j=rng.permutation(i); M[i]=M[j]; Rv[i]=Rv[j]; q[i]=q[j]
    U=rng.normal(size=(n,3)); U/=np.linalg.norm(U,axis=1)[:,None]; d=(rng.uniform(dmin**3,dmax**3,n))**(1/3)
    p0=m.encode(U*d[:,None],M,Rv,q); lam=m.lam if m.lam>0 else 10
    p=relax(m,p0,lam) if (m.L["cont"] or m.L["net"] or m.L["trend"] or not mir) else relax(m,p0,10)
    dec=m.decode(p); X2,M2,Rv2,Rc2,dl2=dec; h,a,s=m.split(dec); w=m.pattern_X(X2,M2)
    return dict(left=float(rms(R-w)/FLOOR),ok=bool(h.max(initial=0)<=1e-3 and a.max(initial=0)<=1e-3),soft_max=float(s.max(initial=0)),
                soft_rms=float(np.sqrt(np.mean(s**2))) if len(s) else 0.0,props=network_props(kinds,X2,M2,Rv2),X=X2.tolist(),M=M2.tolist(),Rv=Rv2.tolist())
if __name__=="__main__":
    npr,K,level,sym,N=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3],sys.argv[4],int(sys.argv[5])
    b=best(level,sym,N); X=np.array(b["X"]); d=np.linalg.norm(X,axis=1); q=[1.0 if kk=="C" else -dl/0.7 for kk,dl in zip(b["kinds"],b["delta"])]
    tasks=[(level,sym=="m",b["kinds"],b["X"],b["M"],b["Rv"],q,max(0.6,d.min()),d.max(),s) for s in range(K)]
    os.makedirs("null",exist_ok=True); t0=time.time()
    with Pool(npr) as pool: out=pool.map(one,tasks,chunksize=2)
    json.dump(dict(fit=dict(left=b["left"],props=b["props"],start=b["start"],soft_rms=b["soft_rms"]),null=out),open(f"null/{level}_{sym}_N{N}.json","w"))
    L=np.array([o["left"] for o in out]); print(level,sym,N,f"{time.time()-t0:.0f}s fit={b['left']:.3f} null min={L.min():.2f} med={np.median(L):.1f} frac<=fit={np.mean(L<=b['left']):.4f} ok={np.mean([o['ok'] for o in out]):.2f} soft_med={np.median([o['soft_rms'] for o in out]):.3f} fit_soft={b['soft_rms']:.3f}")

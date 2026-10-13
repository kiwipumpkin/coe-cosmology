"""usage: python3 run_ep.py NPROC  -> raw/<setup>.json"""
import os,sys,json,time,numpy as np
os.environ.setdefault("OMP_NUM_THREADS","1")
from multiprocessing import Pool
from ep_core import *
sys.path.insert(0,os.path.dirname(ME)); from extend import encode
OFFS=os.environ.get("OUTDIR","raw"); P=os.path.dirname(ME); RES=json.load(open(f"{P}/results.json"))
SETUPS=[]
for sg in ("over","under","free"):
    for n in (1,2,3): SETUPS.append((f"{sg}_{n}p",sg,n,"",24))
for sg in ("over","under","free"):
    SETUPS.append((f"{sg}_2p+4C",sg,2,"CCCC",16)); SETUPS.append((f"{sg}_3p+4C",sg,3,"CCCC",16))
    SETUPS.append((f"{sg}_2p+4C4V",sg,2,"CCCCVVVV",16))
def p0(m,rng,i):
    lo,hi=m.lo,m.hi; p=rng.uniform(lo,hi)
    p[0]=rng.uniform(lo[0],hi[0])
    for k in range(m.np_):
        b=1+5*k; p[b]=rng.uniform(0,0.5); p[b+1]=rng.uniform(-np.pi,np.pi); p[b+2]=np.arcsin(rng.uniform(-1,1)); p[b+3]=rng.uniform(0,1); p[b+4]=rng.uniform(m.pl,m.ph)
    if m.kinds:
        o=1+5*m.np_; n=m.ord.n
        if i%2==0:   # warm: old best for this composition
            br=RES[m.kinds]["best_run"]; p[o:]=encode(m.kinds,br["X"],br["M"],br["Rv"])+rng.normal(0,0.02,4*n)
        else:
            q=np.zeros((n,4)); q[:,0]=rng.uniform(0,0.6,n); q[:,1]=rng.uniform(-np.pi,np.pi,n); q[:,2]=np.arcsin(rng.uniform(-1,1,n)); q[:,3]=rng.uniform(0,1,n); p[o:]=q.ravel()
    return np.clip(p,lo,hi)
def job(a):
    name,sg,n,kinds,i=a; m=EP(sg,n,kinds); rng=np.random.default_rng(abs(hash((name,i,OFFS)))%2**31); t0=time.time()
    p=fit(m,p0(m,rng,i)); p=fit(m,p)
    du,Xp,Mp,Rp,dp,Xo,Mo,Rvo,Rco=m.decode(p); w=m.pattern(p)
    return dict(setup=name,sign=sg,npart=n,kinds=kinds,start=i,left=float(m.leftover(p)),viol=float(m.violations(p).max(initial=0)),
        net=float(m.net(p)),delta_us=float(du),M_sphere=float(du*MS1),Xp=Xp.tolist(),Mp=Mp.tolist(),Rp=Rp.tolist(),dp=dp.tolist(),
        Xo=Xo.tolist(),Mo=Mo.tolist(),Rvo=Rvo.tolist(),p=p.tolist(),amps=(amps(R-w)/FLOOR).tolist(),secs=time.time()-t0)
if __name__=="__main__":
    npr=int(sys.argv[1]); only=sys.argv[2:] ; OUT=os.environ.get("OUTDIR","raw"); os.makedirs(OUT,exist_ok=True)
    tasks=[(nm,sg,n,k,i) for nm,sg,n,k,ns in SETUPS if (not only or nm in only) for i in range(ns)]
    t0=time.time(); res={}
    with Pool(npr) as pool:
        for r in pool.imap_unordered(job,tasks):
            res.setdefault(r["setup"],[]).append(r); ok=r["viol"]<=1e-3 and abs(r["net"])<=NET_TOL
            print(f"{time.time()-t0:6.0f}s {r['setup']:14s} #{r['start']:2d} left={r['left']:.2f} ok={ok} d_us={r['delta_us']:+.3f} net={r['net']:.1e}",flush=True)
            json.dump(res[r["setup"]],open(f"{OUT}/{r['setup']}.json","w"))
    print("done",time.time()-t0)

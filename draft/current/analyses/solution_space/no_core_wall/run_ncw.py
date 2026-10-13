"""usage: NCW=smooth|open python3 run_ncw.py NPROC [asym]  -> raw_<var>/ (or raw_<var>_asym/)"""
import os,sys,json,time,numpy as np
os.environ.setdefault("OMP_NUM_THREADS","1")
from multiprocessing import Pool
from ncw_core import *
P=os.path.dirname(ME); sys.path.insert(0,P); from extend import encode
import fit_run
ASYM=len(sys.argv)>2 and sys.argv[2]=="asym"
if ASYM:
    sys.path.insert(0,os.path.join(P,"asym")); import asym_core as AC
CFG=[("CCCC",8,4),("CCCV",6,6),("CCCCCV",6,6),("CCCCVV",6,6),("CCCCVVVV",8,8),("CCCCCVVV",8,8),("CCCCCCVVVV",8,8),("CCCCCCVVVVV",10,6),("CCCCCCVVVVVV",10,6)]
if ASYM: CFG=[("CCCCCCVVVVV",10,6),("CCCCCCVVVVVV",10,6),("CCCCCCCVVVVV",8,8)]
if ASYM and os.environ.get("ASYM_SMALL"): CFG=[("CCCCCCVVVVV",6,4),("CCCCCCVVVVVV",6,4)]
def sources():
    res=json.load(open(f"{P}/results.json")); ext=json.load(open(f"{P}/extend_mixed.json")); S=[]
    for k,v in res.items():
        if "S" not in k and len(k)>=4: S.append(dict(v["best_run"],kinds=k))
    S+=[s for s in ext]; 
    for f in ("cold_CCCCCCVVVVVV_polished.json","cold_CCCCCCVVVVV_polished.json"): S+=json.load(open(f"{P}/{f}"))
    if ASYM:
        for f in os.listdir(f"{P}/asym/raw"): S+=[r for r in json.load(open(f"{P}/asym/raw/{f}"))["runs"] if r["left"]<5]
    return S
def warm(kinds,S,rng):
    nC,nV=kinds.count("C"),kinds.count("V")
    cand=[s for s in S if s["kinds"].count("C")>=nC and s["kinds"].count("V")>=nV and len(s["kinds"])<=len(kinds)+3] or S
    cand=sorted(cand,key=lambda s:s["left"])[:12]; s=cand[rng.integers(len(cand))]
    k=list(s["kinds"]); X=np.array(s["X"]); M=np.array(s["M"]); Rv=np.array(s["Rv"])
    iC=[i for i in range(len(k)) if k[i]=="C"]; iV=[i for i in range(len(k)) if k[i]=="V"]
    iC=list(rng.permutation(iC))[:nC]; iV=list(rng.permutation(iV))[:nV]
    def fill(ids,n,t):
        Xs=[X[i] for i in ids]; Ms=[M[i] for i in ids]; Rs=[Rv[i] for i in ids]
        while len(Xs)<n:
            d=rng.uniform(16,35); u=rng.normal(size=3); u/=np.linalg.norm(u); Xs.append(d*u); Ms.append(10**rng.uniform(51,52.3) if t=="C" else 0); Rs.append(rng.uniform(1,4))
        return Xs,Ms,Rs
    a=fill(iC,nC,"C"); b=fill(iV,nV,"V")
    Xs=np.array(a[0]+b[0]); Xs=Xs*np.exp(rng.normal(0,0.05,(len(Xs),1)))
    p=encode(kinds,Xs,np.clip(np.array(a[1]+b[1],float),MLO,MHI),np.array(a[2]+b[2]))
    p=p.reshape(-1,4); p[:,1]+=rng.normal(0,0.08,len(p)); return p.ravel()
def cold(n,rng):
    p=np.zeros((n,4)); p[:,0]=rng.uniform(0,0.6,n); p[:,1]=rng.uniform(-np.pi,np.pi,n); p[:,2]=np.arcsin(rng.uniform(-1,1,n)); p[:,3]=rng.uniform(0,1,n); return p.ravel()
def job(a):
    kinds,p0,tag=a; m=NCW(kinds); t0=time.time()
    if ASYM:
        p=AC.fit_one(m,p0); p=AC.fit_one(m,p)
    else:
        p=fit_run.fit_one(m,p0); p=fit_run.fit_one(m,p)
    X,M,Rv,Rc=m.decode(p); w=m.pattern(p)
    out=dict(kinds=kinds,start=tag,var=VAR,left=float(rms(R-w)/FLOOR),viol=float(m.violations(p).max(initial=0)),p=p.tolist(),X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),
             amps=(amps(R-w)/FLOOR).tolist(),core_depth=core_depth(X).tolist(),secs=time.time()-t0)
    if ASYM:
        nm,npl,best,plane=AC.hard_check(m.k,X,M,Rv); out.update(mirror_matches=nm,on_plane=npl,asym_viol=float(AC.asym_viol(m.k,X,M,Rv).max(initial=0)),odd_frac=AC.odd_fraction(m,p))
    return out
if __name__=="__main__":
    npr=int(sys.argv[1]); od=f"raw_{VAR}"+("_asym" if ASYM else ""); os.makedirs(od,exist_ok=True)
    S=sources(); tasks=[]
    for ci,(k,nw,nc) in enumerate(CFG):
        rng=np.random.default_rng(700+ci+(50 if VAR=="open" else 0))
        tasks+=[(k,warm(k,S,rng),f"warm{i}") for i in range(nw)]+[(k,cold(len(k),rng),f"cold{i}") for i in range(nc)]
    tasks.sort(key=lambda t:-len(t[0]))   # long ones first
    res={}; t0=time.time()
    with Pool(npr) as pool:
        for r in pool.imap_unordered(job,tasks):
            res.setdefault(r["kinds"],[]).append(r)
            print(f"{time.time()-t0:6.0f}s {VAR} {r['kinds']:13s} {r['start']:7s} left={r['left']:.2f} viol={r['viol']:.1e} maxdepth={max(r['core_depth']):+.1f}",flush=True)
            json.dump(res[r["kinds"]],open(f"{od}/{r['kinds']}.json","w"))
    print("done",time.time()-t0)

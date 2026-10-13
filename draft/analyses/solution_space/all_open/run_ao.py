"""usage: AO=smooth|open python3 run_ao.py NPROC NSTARTS [kinds...]  -> raw_<var>/KINDS.json   (env ASYM=1 -> no-mirror rule, raw_<var>_asym/)"""
import os,sys,json,time,glob,numpy as np
os.environ.setdefault("OMP_NUM_THREADS","1")
from multiprocessing import Pool
from ao_core import *
ASYM=bool(os.environ.get("ASYM"))
if ASYM:
    sys.path.insert(0,os.path.join(SS,"asym")); import asym_core as AC
sys.path.insert(0,SS); import fit_run
CFG=["CCCC","CCCV","CCCCCC","CCCCCV","CCCCVV","CCCCCCCC","CCCCVVVV","CCCCCVVV","CCCCCCCCCC","CCCCCCVVVV","CCCCCCCCCCC","CCCCCCVVVVV","CCCCCCCCCCCC","CCCCCCVVVVVV"]
def sources(kinds):
    S=[]
    f=f"{CS}/inside_visible/raw_cap2e52/{kinds}.json"
    if os.path.exists(f): S+=[dict(r,src="inside") for r in sorted([r for r in json.load(open(f))["runs"] if r["viol"]<=1e-3],key=lambda r:r["left"])[:4]]
    f=f"{SS}/no_core_wall/raw_smooth/{kinds}.json"
    if os.path.exists(f): S+=[dict(r,src="outside") for r in sorted([r for r in json.load(open(f)) if r["viol"]<=1e-3],key=lambda r:r["left"])[:4]]
    return S
def cold(n,rng):
    p=np.zeros((n,4)); p[:,0]=rng.uniform(0,1,n)**1.5; p[:,1]=rng.uniform(-np.pi,np.pi,n); p[:,2]=np.arcsin(rng.uniform(-1,1,n)); p[:,3]=rng.uniform(0,1,n); return p.ravel()
def mixed(kinds,rng):
    """half the objects from the best inside fit, half from the best outside fit (same kinds) -- if both exist"""
    S=sources(kinds); I=[s for s in S if s["src"]=="inside"]; O=[s for s in S if s["src"]=="outside"]
    if not I or not O: return None
    a,b=I[0],O[0]; k=list(kinds); X=[];M=[];Rv=[]
    for i in range(len(k)):
        s=a if rng.uniform()<0.5 else b; X.append(s["X"][i]); M.append(s["M"][i]); Rv.append(s["Rv"][i])
    return encode(kinds,X,M,Rv)
def job(a):
    kinds,p0,tag=a; m=AO(kinds); t0=time.time()
    if ASYM: p=AC.fit_one(m,p0); p=AC.fit_one(m,p)
    else: p=fit_run.fit_one(m,p0); p=fit_run.fit_one(m,p)
    X,M,Rv,Rc=m.decode(p); w=m.pattern(p)
    out=dict(kinds=kinds,start=tag,var=VAR,left=float(rms(R-w)/FLOOR),viol=float(m.violations(p).max(initial=0)),p=p.tolist(),X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),Rc=Rc.tolist(),
             amps=(amps(R-w)/FLOOR).tolist(),secs=time.time()-t0)
    if ASYM:
        nm,npl,best,plane=AC.hard_check(m.k,X,M,Rv); out.update(mirror_matches=nm,on_plane=npl,odd_frac=AC.odd_fraction(m,p))
    return out
if __name__=="__main__":
    npr,ns=int(sys.argv[1]),int(sys.argv[2]); cfg=sys.argv[3:] or CFG
    od=f"raw_{VAR}"+("_asym" if ASYM else ""); os.makedirs(od,exist_ok=True); tasks=[]
    for ci,k in enumerate(cfg):
        rng=np.random.default_rng(900+ci+(77 if VAR=="open" else 0)+(500 if ASYM else 0)); S=sources(k); n=len(k)
        for i in range(ns):
            c=i%4
            if c==0 and S: s=S[(i//4)%len(S)]; p=encode(k,s["X"],s["M"],s["Rv"]).reshape(-1,4); p[:,1]+=rng.normal(0,0.1 if ASYM else 0.03,n); p=p.ravel(); tag=f"warm_{s['src']}{i}"
            elif c==1 and mixed(k,rng) is not None: p=mixed(k,rng); tag=f"mixed{i}"
            else: p=cold(n,rng); tag=f"cold{i}"
            tasks.append((k,p,tag))
    tasks.sort(key=lambda t:-len(t[0])); res={}; t0=time.time()
    with Pool(npr) as pool:
        for r in pool.imap_unordered(job,tasks):
            res.setdefault(r["kinds"],[]).append(r); d=np.linalg.norm(r["X"],axis=1); nin=int((d<RV).sum())
            print(f"{time.time()-t0:6.0f}s {VAR} {r['kinds']:13s} {r['start']:14s} left={r['left']:.3f} viol={r['viol']:.1e} inside={nin}/{len(d)}",flush=True)
            json.dump(res[r["kinds"]],open(f"{od}/{r['kinds']}.json","w"))
    print("done",time.time()-t0)

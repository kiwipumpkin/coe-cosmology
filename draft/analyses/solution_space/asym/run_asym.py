"""usage: python3 run_asym.py KINDS NWARM NCOLD SEED [NPROC] -> raw/KINDS.json (env ASYM_ANG=10 -> raw10/)"""
import os, sys, json, time, numpy as np
os.environ.setdefault("OMP_NUM_THREADS","1")
from multiprocessing import Pool
from asym_core import *
ME=os.path.dirname(os.path.abspath(__file__)); P=os.path.dirname(ME); sys.path.insert(0,P); from extend import encode

def load_sources():
    ext=json.load(open(f"{P}/extend_mixed.json")); res=json.load(open(f"{P}/results.json")); src=[]
    src.append(ext[-1])                       # old floor solution (6C+5V, mirror-symmetric)
    src+=ext[1:-1]
    for k in ("CCCCVVVV","CCCCCVVV","CCCCCCVV","CCCVVVVV"): src.append(dict(res[k]["best_run"],kinds=k))
    for f in ("cold_CCCCCCVVVVVV_polished.json","cold_CCCCCCVVVVV_polished.json"):
        src+=sorted(json.load(open(f"{P}/{f}")),key=lambda r:r["left"])[:3]
    return src
def obj_list(s):
    k=list(s["kinds"]); return [(k[i],np.array(s["X"][i]),s["M"][i],s["Rv"][i]) for i in range(len(k))]
def rand_obj(t,rng):
    d=rng.uniform(16,35); l=rng.uniform(0,2*np.pi); b=np.arcsin(rng.uniform(-1,1))
    X=d*np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
    return (t,X,10**rng.uniform(51,52.3),rng.uniform(1,4))
def warm_p0(kinds,src,rng):
    nC,nV=kinds.count("C"),kinds.count("V"); ob=obj_list(src)
    Cs=[o for o in ob if o[0]=="C"]; Vs=[o for o in ob if o[0]=="V"]
    def pick(L,n,t):
        L=[L[i] for i in rng.permutation(len(L))][:n]
        while len(L)<n: L.append(rand_obj(t,rng))
        return L
    objs=pick(Cs,nC,"C")+pick(Vs,nV,"V"); X=[];M=[];Rv=[]
    for t,x,m,rv in objs:
        d=np.linalg.norm(x); l=np.arctan2(x[1],x[0]); b=np.arcsin(x[2]/d)
        l+=np.radians(rng.normal(0,15)); b=np.clip(b+np.radians(rng.normal(0,8)),-1.4,1.4); d*=np.exp(rng.normal(0,0.10))
        # push off the mirror plane
        u=np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
        if np.degrees(np.arcsin(abs(u@NRM)))<PLANE_TOL*1.5: l+=np.radians(rng.choice([-1,1])*rng.uniform(1.5*PLANE_TOL,3*PLANE_TOL))
        X.append(d*np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)]))
        M.append(m*np.exp(rng.normal(0,0.3)) if t=="C" else m); Rv.append(np.clip(rv*np.exp(rng.normal(0,0.15)),0.5,4) if t=="V" else 0)
    return encode(kinds,X,np.clip(np.array(M,float),MLO,MHI) if True else M,Rv)
def cold_p0(n,rng):
    p=np.zeros((n,4)); p[:,0]=rng.uniform(0,0.6,n); p[:,1]=rng.uniform(-np.pi,np.pi,n); p[:,2]=np.arcsin(rng.uniform(-1,1,n)); p[:,3]=rng.uniform(0,1,n); return p.ravel()
def job(a):
    kinds,p0,tag=a; m=Model(kinds); t0=time.time(); p=fit_one(m,p0)
    # polish: 2 more passes from the result
    p=fit_one(m,p)
    X,M,Rv,Rc=m.decode(p); w=m.pattern(p); k=m.k
    nm,npl,best,plane=hard_check(k,X,M,Rv)
    return dict(kinds=kinds,start=tag,left=float(rms(R-w)/FLOOR),viol=float(m.violations(p).max(initial=0)),
        asym_viol=float(asym_viol(k,X,M,Rv).max(initial=0)),mirror_matches=nm,on_plane=npl,
        closest_pair=[float(v) for v in best[2:]] if best else None,min_plane_deg=float(plane.min()),odd_frac=odd_fraction(m,p),
        p=p.tolist(),X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),amps=(amps(R-w)/FLOOR).tolist(),secs=time.time()-t0)
if __name__=="__main__":
    kinds,nw,nc,seed=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]); npr=int(sys.argv[5]) if len(sys.argv)>5 else 4
    out_dir="raw" if ANG_TOL==20 else f"raw{int(ANG_TOL)}"; os.makedirs(out_dir,exist_ok=True)
    rng=np.random.default_rng(seed); src=load_sources(); tasks=[]
    for i in range(nw):
        s=src[0] if i%2==0 else src[1+rng.integers(len(src)-1)]
        tasks.append((kinds,warm_p0(kinds,s,rng),f"warm{i}"))
    for i in range(nc): tasks.append((kinds,cold_p0(len(kinds),rng),f"cold{i}"))
    t0=time.time()
    with Pool(npr) as pool: out=pool.map(job,tasks,chunksize=1)
    json.dump(dict(kinds=kinds,ang_tol=ANG_TOL,runs=out,secs=time.time()-t0),open(f"{out_dir}/{kinds}.json","w"))
    ok=[r for r in out if r["viol"]<=1e-3 and r["mirror_matches"]==0 and r["on_plane"]==0]
    L=sorted(r["left"] for r in ok)
    print(kinds,f"{time.time()-t0:.0f}s ok={len(ok)}/{len(out)} best={L[0] if L else None} top5={np.round(L[:5],2).tolist()}",flush=True)

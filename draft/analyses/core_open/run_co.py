"""usage: python3 run_co.py NPROC [scenario ...]   scenarios: out asym ao ep  -> raw/<scen>_<kinds>.json"""
import os,sys,json,time,glob,numpy as np
os.environ.setdefault("OMP_NUM_THREADS","1")
from multiprocessing import Pool
from co_core import *
import fit_run
sys.path.insert(0,os.path.join(SS,"asym")); import asym_core as AC
OLD=SS
CFG={"out":[("CCCC",12),("CCCCVV",12),("CCCCVVVV",12),("CCCCCVVVV",10),("CCCCCCVVVV",10),("CCCCCCVVVVV",10),("CCCCCCVVVVVV",10)],
     "asym":[("CCCCCCVVVV",8),("CCCCCCVVVVV",8),("CCCCCCVVVVVV",8)],
     "ao":[("CCCC",8),("CCCCVV",8),("CCCCCC",8),("CCCCCCC",8),("CCCCVVVV",8),("CCCCCCCC",8)],
     "ep":[("free_1p",8),("free_2p",8),("free_3p",8),("free_2p+CCCC",8),("free_2p+CCCCVVVV",10),("over_2p+CCCCVVVV",8),("free_2p+CCCCCVVVV",8),("free_2p+CCCCCCVVVV",8)]}
def _ld(f):
    if not os.path.exists(f): return []
    r=json.load(open(f)); r=r["runs"] if isinstance(r,dict) else r
    return sorted([x for x in r if x.get("viol",0)<=1e-3],key=lambda x:x["left"])
def sources(scen,kinds):
    if scen=="out": return _ld(f"{OLD}/no_core_wall/raw_open/{kinds}.json")[:3]+_ld(f"{OLD}/raw/{kinds}.json")[:2]
    if scen=="asym": return _ld(f"{OLD}/no_core_wall/raw_smooth_asym/{kinds}.json")[:3]+_ld(f"{OLD}/asym/raw/{kinds}.json")[:3]
    if scen=="ao": return _ld(f"{OLD}/all_open/raw_open/{kinds}.json")[:4]
    return []
def model(scen,name):
    if scen=="ep":
        sg,rest=name.split("_"); npart=int(rest[0]); kinds=rest.split("+")[1] if "+" in rest else ""
        return COEP(sg,npart,kinds)
    return COAO(name) if scen=="ao" else CO(name)
def body_seed(m,rng,k,in_core_flag):
    """one random body param row; in core if flag"""
    if k=="C": M=10**rng.uniform(51,np.log10(MHI)); Rc=(3*M/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3); Rv=1.0
    else: Rv=rng.uniform(1,4); Rc=Rv; M=-1
    if in_core_flag: X=random_core_point(rng,Rc)
    else:
        v=rng.normal(size=3); v/=np.linalg.norm(v); X=v*rng.uniform(16,35)
    return X,M,Rv
def start(scen,name,m,rng,i,S):
    """tags: warm (old best), cold (uniform in the open reach), core (part or all bodies moved into the core)"""
    if scen=="ep":
        p=rng.uniform(m.lo,m.hi)
        for q in range(m.np_):
            b=1+5*q; p[b]=rng.uniform(0,1); p[b+1]=rng.uniform(-np.pi,np.pi); p[b+2]=np.arcsin(rng.uniform(-1,1)); p[b+3]=rng.uniform(0,1); p[b+4]=rng.uniform(m.pl,m.ph)
        mode=["warm","cold","core"][i%3]; p[0]=rng.uniform(-0.05,0.05) if mode!="cold" else p[0]
        if mode=="core":     # partners seeded in the core: centre in the CoE interior
            for q in range(m.np_):
                if q==0 or rng.uniform()<0.5:
                    b=1+5*q; X=random_core_point(rng,5.0); d=np.linalg.norm(X); U=X/d; Rp=1+p[b+3]*14
                    dmn=RMIN+Rp; dmx=dmax_dir(U[None],np.array([Rp]),60.0)[0]
                    p[b]=np.clip((d-dmn)/(dmx-dmn),0,1); p[b+1]=np.arctan2(U[1],U[0]); p[b+2]=np.arcsin(U[2])
        if m.kinds:
            o=1+5*m.np_; om=m.ord; k=list(m.kinds)
            base=None
            if mode=="warm":
                prev=json.load(open(f"{OLD}/results.json")); kk=m.kinds
                if kk in prev: br=prev[kk]["best_run"]; base=(br["X"],br["M"],br["Rv"])
                else:
                    ext=[e for e in json.load(open(f"{OLD}/extend_mixed.json")) if e["kinds"]==kk]
                    if ext: br=min(ext,key=lambda e:e["left"]); base=(br["X"],br["M"],br["Rv"])
            Xs=[];Ms=[];Rs=[]
            for j,t in enumerate(k):
                if base is not None: Xs.append(base[0][j]); Ms.append(abs(base[1][j]) if t=="C" else 1); Rs.append(base[2][j])
                else:
                    X,M,Rv=body_seed(om,rng,t,mode=="core" and rng.uniform()<0.5); Xs.append(X); Ms.append(abs(M)); Rs.append(Rv)
            p[o:]=om.encode(Xs,Ms,Rs)
        return np.clip(p,m.lo,m.hi),mode
    k=list(name); n=len(k); c=i%4
    if c==0 and S:
        s=S[(i//4)%len(S)]; p=m.encode(s["X"],s["M"],s["Rv"]).reshape(-1,4); p[:,1]+=rng.normal(0,0 if i==0 else (0.01 if scen=="ao" else 0.03),n); return p.ravel(),"warm"
    if c==1:
        p=np.zeros((n,4)); p[:,0]=rng.uniform(0,1,n); p[:,1]=rng.uniform(-np.pi,np.pi,n); p[:,2]=np.arcsin(rng.uniform(-1,1,n)); p[:,3]=rng.uniform(0,1,n); return p.ravel(),"cold"
    # core-seeded: fraction 1/3, 1/2 or all of the bodies placed at random points in the core, rest from a warm start
    frac=[1/3,1/2,1.0][(i//4)%3] if c==2 else [1.0,1/2,1/3][(i//4)%3]
    s=S[rng.integers(len(S))] if S else None; Xs=[];Ms=[];Rs=[]
    ncore=max(1,int(round(frac*n))); idx=set(rng.permutation(n)[:ncore].tolist())
    for j,t in enumerate(k):
        if j in idx or s is None: X,M,Rv=body_seed(m,rng,t,j in idx); Xs.append(X); Ms.append(abs(M)); Rs.append(Rv)
        else: Xs.append(s["X"][j]); Ms.append(abs(s["M"][j]) if t=="C" else 1); Rs.append(s["Rv"][j])
    return m.encode(Xs,Ms,Rs),f"core{ncore}of{n}"
def job(a):
    scen,name,p0,tag=a; m=model(scen,name); t0=time.time()
    if scen=="ep":
        p=ep_core.fit(m,p0); p=ep_core.fit(m,p)
        du,Xp,Mp,Rp,dp,Xo,Mo,Rvo,Rco=m.decode(p); w=m.pattern(p)
        X=np.r_[Xp,Xo]; out=dict(delta_us=float(du),net=float(m.net(p)),Xp=Xp.tolist(),Mp=Mp.tolist(),Rp=Rp.tolist(),dp=dp.tolist(),
            Xo=Xo.tolist(),Mo=Mo.tolist(),Rvo=Rvo.tolist())
        ok=m.violations(p).max(initial=0)<=1e-3 and abs(m.net(p))<=ep_core.NET_TOL
    else:
        if scen=="asym": p=AC.fit_one(m,p0); p=AC.fit_one(m,p)
        else: p=fit_run.fit_one(m,p0); p=fit_run.fit_one(m,p)
        X,M,Rv,Rc=m.decode(p); w=m.pattern(p); out=dict(X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),Rc=Rc.tolist())
        ok=m.violations(p).max(initial=0)<=1e-3
        if scen=="asym":
            nm,npl,best,plane=AC.hard_check(m.k,X,M,Rv); av=float(AC.asym_viol(m.k,X,M,Rv).max(initial=0))
            out.update(mirror_matches=nm,on_plane=npl,asym_viol=av); ok=ok and nm==0 and npl==0
    left=float(rms(R-w)/FLOOR)
    out.update(scen=scen,name=name,start=tag,left=left,ok=bool(ok),viol=float(m.violations(p).max(initial=0)),p=p.tolist(),
               resid=(R-w).tolist(),amps=(amps(R-w)/FLOOR).tolist(),core_depth=core_depth(X).tolist(),dist=np.linalg.norm(X,axis=1).tolist(),secs=time.time()-t0)
    return out
if __name__=="__main__":
    npr=int(sys.argv[1]); scens=sys.argv[2:] or ["out","asym","ao","ep"]; os.makedirs("raw",exist_ok=True); tasks=[]
    for si,sc in enumerate(scens):
        for ci,(name,ns) in enumerate(CFG[sc]):
            rng=np.random.default_rng(4242+100*si+ci); m=model(sc,name); S=sources(sc,name) if sc!="ep" else []
            for i in range(ns):
                p0,tag=start(sc,name,m,rng,i,S); tasks.append((sc,name,p0,f"{tag}{i}"))
    def cost(t):
        m=model(t[0],t[1]); return -len(m.lo)*(3 if t[0]=="asym" else 1)
    tasks.sort(key=cost); res={}; t0=time.time()
    print("tasks",len(tasks),flush=True)
    with Pool(npr) as pool:
        for r in pool.imap_unordered(job,tasks):
            key=f"{r['scen']}_{r['name']}"; res.setdefault(key,[]).append(r); ncore=int((np.array(r["core_depth"])>0).sum())
            print(f"{time.time()-t0:6.0f}s {key:28s} {r['start']:12s} left={r['left']:.3f} ok={r['ok']} core={ncore}/{len(r['core_depth'])} maxdepth={max(r['core_depth']):+.1f} t={r['secs']:.0f}",flush=True)
            json.dump(res[key],open(f"raw/{key}.json","w"))
    print("done",time.time()-t0)

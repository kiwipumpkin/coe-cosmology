"""usage: python3 run_shallow.py near <nseeds> | swap | outer <random|foam> <nreal> | grow <foam|random|cmb_aware> <seed>"""
import sys, os, json, time, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from shallow_core import *
HERE=os.path.dirname(os.path.abspath(__file__)); SEQ=os.path.dirname(HERE)
RADII=(1.5,2.25,3.0,3.75)

def near(seed,nbank=5000,max_voids=48):
    rng=np.random.default_rng(seed); pk=Packer(h=4.0); w=np.zeros(360); traj=[float(ldiff(w))]; placed=[]
    banks={}
    for Rv in RADII:
        P=rand_in_shell(rng,RMIN,LOCAL_RMAX,nbank); P=P[region_ok(P,Rv,rmax=LOCAL_RMAX)]
        banks[Rv]=(P,bank_patterns(P,Rv))
    for k in range(max_voids):
        res=R-w; base=res@res; best=None
        for Rv,(P,F) in banks.items():
            S=((res[None]-F)**2).sum(1)
            for i in np.argsort(S)[:300]:
                if S[i]>=base*(1-1e-5): break
                if pk.free(P[i],Rv):
                    if best is None or S[i]<best[0]: best=(S[i],P[i],Rv,F[i])
                    break
        if best is None: break
        _,x,Rv,f=best; pk.add(x,Rv,-float(def_mass(x,Rv)[0]),"void"); w=w+f; traj.append(float(ldiff(w))); placed.append((x.tolist(),Rv))
    X,r,M,_=pk.arrays()
    return dict(seed=seed,left=traj[-1],traj=traj,X=X.tolist(),Rv=r.tolist(),n=len(r),min_gap=float(check_no_overlap(X,r,0.0)) if len(r)>1 else None,
                amps=(amps(R-w)/FLOOR)[:6].tolist(),w=w.tolist())

def swap():
    z=np.load(os.path.join(SEQ,"local_frozen.npz"),allow_pickle=True); X,r,M,k=z["X"],z["r"],z["M"],z["kind"].astype(str)
    v=k=="void"
    w_s=shallow_pattern(X[v],r[v])                       # voids -> shallow profile (self-compensated), clusters removed
    w_c=pattern(X[~v],M[~v])                             # old compact clusters, unchanged
    w_old=pattern(X,M)
    return dict(deep_original=float(ldiff(w_old)),shallow_voids_no_clusters=float(ldiff(w_s)),
                shallow_voids_plus_old_clusters=float(ldiff(w_s+w_c)),n_voids=int(v.sum()),n_clusters=int((~v).sum()))

def fill(cfg,rng,Rv=1.5,r0=LOCAL_RMAX,r1=S_MAX_G,pk=None,max_fail=3000,max_voids=10**9,rule_cb=None):
    """sequential void cores (no overlap), random drop or grown next to neighbours (foam)."""
    pk=pk or Packer(h=4.0); new=[]; fronts=[]; fail=0
    while fail<max_fail and len(new)<max_voids:
        if cfg=="foam" and fronts:
            j=rng.integers(len(fronts)); p=fronts[j]; ok=False
            for _ in range(24):
                x=p+rand_dir(rng)*2*Rv*rng.uniform(1.02,1.15)
                if region_ok(x,Rv,rmin=r0,rmax=r1)[0] and pk.free(x,Rv): ok=True; break
            if not ok: fronts.pop(j); continue
        else:
            x=rand_in_shell(rng,r0,r1,1)[0]
            if not (region_ok(x,Rv,rmin=r0,rmax=r1)[0] and pk.free(x,Rv)): fail+=1; continue
        pk.add(x,Rv,0.0,"void"); new.append(x); fronts.append(x); fail=0 if cfg=="foam" else fail
    return np.array(new),pk

def outer(cfg,nreal):
    nf=json.load(open(os.path.join(HERE,"near_fit.json"))); b=nf["best"]
    Xn=np.array(b["X"]); rn=np.array(b["Rv"]); w0=np.array(b["w"]); L0=float(ldiff(w0)); rows=[]
    for i in range(nreal):
        rng=np.random.default_rng([77,i]); pk=Packer(h=4.0)
        for x,r in zip(Xn,rn): pk.add(x,r,0.0,"void")
        Xo,pk=fill(cfg,rng,pk=pk)
        wo=shallow_pattern(Xo,1.5) if len(Xo) else np.zeros(360)
        X,r,_,_=pk.arrays()
        rows.append(dict(units=len(Xo),left=float(ldiff(w0+wo)),outer_rms=float(rms(wo)/FLOOR),min_gap=float(check_no_overlap(X,r,0.0))))
        print(cfg,i,rows[-1],flush=True)
    json.dump(dict(cfg=cfg,frozen=L0,rows=rows),open(os.path.join(HERE,f"outer_{cfg}.json"),"w"),indent=1)

def grow(rule,seed,Rv=1.5):
    rng=np.random.default_rng(seed); pk=Packer(h=4.0); fronts=[]; order=[]; w=np.zeros(360); fail=0; t0=time.time()
    while fail<3 and len(order)<6000:
        cands=[]
        for _ in range(600):
            if rule!="random" and fronts:
                j=rng.integers(len(fronts)); x=fronts[j]+rand_dir(rng)*2*Rv*rng.uniform(1.02,1.12)
            else: x=rand_in_shell(rng,RMIN,S_MAX_G,1)[0]
            if region_ok(x,Rv)[0] and pk.free(x,Rv): cands.append(x)
            if len(cands)>=(30 if rule=="cmb_aware" else 1): break
        if not cands:
            fronts=fronts[len(fronts)//2:]; fail+=1; continue
        fail=0
        if rule=="cmb_aware":
            F=bank_patterns(np.array(cands),Rv); i=int(np.argmin(((R-w)[None]-F)**2).sum(1)) if False else int(np.argmin((((R-w)[None]-F)**2).sum(1))); x=cands[i]; w=w+F[i]
        else: x=cands[0]
        pk.add(x,Rv,0.0,"void"); fronts.append(x); order.append(x)
        if len(order)%1000==0: print(rule,seed,len(order),f"{time.time()-t0:.0f}s",flush=True)
    X=np.array(order); left=np.empty(len(X)); w=np.zeros(360)
    F=bank_patterns(X,Rv)
    C=np.cumsum(F,0); left=np.sqrt(np.mean((R[None]-C)**2,1))/FLOOR
    Xa,ra,_,_=pk.arrays()
    out=dict(rule=rule,seed=seed,n_void=len(X),void_vol_frac=float(len(X)*vol(Rv)/268960.0),final=float(left[-1]),min=float(left.min()),
             start=float(BASE0),min_gap=float(check_no_overlap(Xa,ra,0.0)),traj=left[::5].round(3).tolist(),
             n_near_ring=int(np.sum(ring_d(X)<3*Rv)))
    json.dump(out,open(os.path.join(HERE,f"grow_{rule}_s{seed}.json"),"w"))
    print("DONE",{k:v for k,v in out.items() if k!="traj"},flush=True)

def ring_d(X): return np.sqrt((np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2)

if __name__=="__main__":
    m=sys.argv[1]
    if m=="near":
        runs=[]
        for s in range(int(sys.argv[2])):
            r=near(s); runs.append(r); print("seed",s,"left",round(r["left"],2),"voids",r["n"],"gap",r["min_gap"],flush=True)
        b=min(runs,key=lambda r:r["left"])
        json.dump(dict(best=b,runs=[{k:v for k,v in r.items() if k not in("w",)} for r in runs]),open(os.path.join(HERE,"near_fit.json"),"w"))
    elif m=="swap": d=swap(); print(d); json.dump(d,open(os.path.join(HERE,"swap.json"),"w"),indent=1)
    elif m=="outer": outer(sys.argv[2],int(sys.argv[3]))
    elif m=="grow": grow(sys.argv[2],int(sys.argv[3]))

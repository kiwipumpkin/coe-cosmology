"""Keep adding objects beyond N=8 (warm starts) until leftover < 1x floor, N=20, or 3 stalled steps.
usage: python3 extend.py TRACK   (TRACK = mixed | masses)  -> extend_TRACK.json"""
import os, sys, json, time, numpy as np
os.environ.setdefault("OMP_NUM_THREADS","1")
from multiprocessing import Pool
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from ss_core import *
from fit_run import fit_one
def encode(kinds,X,M,Rv):
    m=Model(kinds); X=np.array(X); M=np.array(M); Rv=np.array(Rv); k=np.array(list(kinds)); p=np.zeros((len(k),4))
    d=np.linalg.norm(X,axis=1); p[:,1]=np.arctan2(X[:,1],X[:,0]); p[:,2]=np.arcsin(X[:,2]/d)
    isC=k=="C"
    p[isC,3]=(np.log10(M[isC])-np.log10(MLO))/np.log10(MHI/MLO); p[~isC,3]=(Rv[~isC]-RLO)/(RHI-RLO)
    p[:,3]=np.clip(p[:,3],0,1)
    Rc=np.where(isC,(3*np.abs(M)/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3),Rv)
    p[:,0]=np.clip((d-RMIN-Rc)/(RMAX-RMIN-Rc),0,1)
    return p.ravel()
def job(a):
    kinds,p0,seed=a; m=Model(kinds); p=fit_one(m,p0); X,M,Rv,Rc=m.decode(p); w=m.pattern(p)
    return dict(kinds=kinds,seed=seed,left=float(rms(R-w)/FLOOR),viol=float(m.violations(p).max(initial=0)),p=p.tolist(),
                X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),amps=(amps(R-w)/FLOOR).tolist())
def new_obj(rng):
    return [rng.uniform(0,0.6),rng.uniform(-np.pi,np.pi),np.arcsin(rng.uniform(-1,1)),rng.uniform(0,1)]
if __name__=="__main__":
    track=sys.argv[1]; K=int(sys.argv[2]) if len(sys.argv)>2 else 4
    start={"mixed":"CCCCVVVV","masses":"CCCCCCCC"}[track]
    b=json.load(open(f"raw/{start}.json")); b=min([r for r in b["runs"] if r["viol"]<=1e-3],key=lambda r:r["left"])
    kinds=start; p=encode(kinds,b["X"],b["M"],b["Rv"]); m=Model(kinds)
    print(track,"start",kinds,"stored",round(b["left"],2),"re-encoded",round(m.leftover(p),2),flush=True)
    steps=[dict(N=len(kinds),kinds=kinds,left=m.leftover(p),amps=b["amps"],X=b["X"],M=b["M"],Rv=b["Rv"],p=p.tolist())]
    rng=np.random.default_rng(7 if track=="mixed" else 8); stall=0; t0=time.time()
    types="CV" if track=="mixed" else "C"
    with Pool(8) as pool:
        while True:
            prev=steps[-1]; N=prev["N"]+1
            if N>20: break
            tasks=[]
            for t in types:
                for k in range(K if track=="mixed" else 2*K):
                    tasks.append((prev["kinds"]+t,np.r_[np.array(prev["p"]),new_obj(rng)],k))
            out=[r for r in pool.map(job,tasks) if r["viol"]<=1e-3]
            best=min(out,key=lambda r:r["left"])
            byt={t:min([r["left"] for r in out if r["kinds"][-1]==t],default=None) for t in types}
            k_=np.array(list(best["kinds"])); M=np.array(best["M"]); Rv=np.array(best["Rv"])
            capM=int(((k_=="C")&(M>=2e52*0.98)).sum()); capV=int(((k_=="V")&(Rv>=3.98)).sum())
            rel=(prev["left"]-best["left"])/prev["left"]
            st=dict(N=N,kinds=best["kinds"],left=best["left"],new_mass_best=byt.get("C"),new_void_best=byt.get("V"),
                    amps=best["amps"],X=best["X"],M=best["M"],Rv=best["Rv"],p=best["p"],caps_mass=capM,caps_void=capV,
                    M_clusters=float(M[k_=="C"].sum()),M_voids=float(M[k_=="V"].sum()),secs=time.time()-t0,n_tasks=len(tasks),n_ok=len(out))
            if best["left"]>=prev["left"]:   # no improvement: keep previous layout + inert-ish? record and count stall
                st.update(note="no improvement over N-1")
            steps.append(st)
            stall=stall+1 if rel<0.02 else 0
            top=np.argsort(best["amps"])[::-1][:3]+1
            print(f"N={N} {best['kinds']} left={best['left']:.2f}x (mass-add {byt.get('C')}, void-add {byt.get('V')}) top modes m={top.tolist()} "
                  f"caps M:{capM} V:{capV} t={time.time()-t0:.0f}s",flush=True)
            json.dump(steps,open(f"extend_{track}.json","w"))
            if best["left"]<1.0 or stall>=3: break

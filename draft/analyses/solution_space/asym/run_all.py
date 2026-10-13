import os,sys,json,time,numpy as np
from multiprocessing import Pool
from run_asym import *
CFG=[("CCCCCCVVVVV",14,10),("CCCCCCVVVVVV",14,10),("CCCCVVVV",10,6),("CCCCCVVV",10,6),("CCCCCVVVV",10,6),("CCCCCCVVV",10,6),
     ("CCCCCCVVVV",10,6),("CCCCCVVVVV",10,6),("CCCCCVVVVVV",10,6),("CCCCCCCVVVVV",10,6)]
if __name__=="__main__":
    npr=int(sys.argv[1]); os.makedirs("raw",exist_ok=True); src=load_sources(); tasks=[]
    for ci,(k,nw,nc) in enumerate(CFG):
        rng=np.random.default_rng(100+ci)
        for i in range(nw):
            s=src[0] if i%2==0 else src[1+rng.integers(len(src)-1)]
            tasks.append((k,warm_p0(k,s,rng),f"warm{i}"))
        for i in range(nc): tasks.append((k,cold_p0(len(k),rng),f"cold{i}"))
    t0=time.time(); res={}
    with Pool(npr) as pool:
        for r in pool.imap_unordered(job,tasks):
            res.setdefault(r["kinds"],[]).append(r)
            ok=r["viol"]<=1e-3 and r["mirror_matches"]==0 and r["on_plane"]==0
            print(f"{time.time()-t0:6.0f}s {r['kinds']:13s} {r['start']:7s} left={r['left']:.2f} ok={ok} odd={r['odd_frac']:.3f}",flush=True)
            for k,v in res.items(): json.dump(dict(kinds=k,ang_tol=ANG_TOL,runs=v),open(f"raw/{k}.json","w"))
    print("done",time.time()-t0)

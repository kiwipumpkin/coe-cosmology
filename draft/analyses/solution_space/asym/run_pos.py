import os,sys
os.environ["ASYM_MODE"]="pos"
import json,time,numpy as np
from multiprocessing import Pool
import run_asym as ra
CFG=[("CCCCCCVVVVV",10,6),("CCCCCCVVVVVV",10,6)]
if __name__=="__main__":
    os.makedirs("raw_pos",exist_ok=True); src=ra.load_sources(); tasks=[]
    for ci,(k,nw,nc) in enumerate(CFG):
        rng=np.random.default_rng(500+ci)
        for i in range(nw): tasks.append((k,ra.warm_p0(k,src[0] if i%2==0 else src[1+rng.integers(len(src)-1)],rng),f"warm{i}"))
        for i in range(nc): tasks.append((k,ra.cold_p0(len(k),rng),f"cold{i}"))
    res={}
    with Pool(int(sys.argv[1])) as pool:
        for r in pool.imap_unordered(ra.job,tasks):
            res.setdefault(r["kinds"],[]).append(r); print(r["kinds"],r["start"],round(r["left"],2),r["mirror_matches"],r["on_plane"],round(r["odd_frac"],3),flush=True)
            json.dump(dict(kinds=r["kinds"],mode="pos",runs=res[r["kinds"]]),open(f"raw_pos/{r['kinds']}.json","w"))

import os,sys,json,numpy as np
from multiprocessing import Pool
sys.path.insert(0,"."); from extend import job
def rep(r):
    for _ in range(4): r=job((r["kinds"],np.array(r["p"]),r["seed"]))
    return r
if __name__=="__main__":
    f=sys.argv[1]; d=json.load(open(f)); good=[r for r in d if r["left"]<5]
    with Pool(8) as pool: out=pool.map(rep,good)
    json.dump(out,open(f.replace(".json","_polished.json"),"w"))
    print(sorted(round(r["left"],3) for r in out))

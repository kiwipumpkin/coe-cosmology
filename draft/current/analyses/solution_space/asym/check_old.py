import json,sys,numpy as np
from asym_core import *
sys.path.insert(0,".."); from extend import encode
ext=json.load(open("../extend_mixed.json")); res=json.load(open("../results.json"))
cands=[("N=%d ext"%s["N"],s["kinds"],s) for s in ext]
for k in ["CCCCVVVV","CCCCCVVV","CCCCVVV"]: cands.append((k,k,res[k]["best_run"]))
for f in ["../cold_CCCCCCVVVVVV_polished.json","../cold_CCCCCCVVVVV_polished.json"]:
    b=min(json.load(open(f)),key=lambda r:r["left"]); cands.append((f.split("_")[2],b["kinds"],b))
for name,k,s in cands:
    kk=np.array(list(k)); X=np.array(s["X"]); M=np.array(s["M"]); Rv=np.array(s["Rv"])
    m=Model(k); p=encode(k,X,M,Rv)
    nm,npl,best,pl=hard_check(kk,X,M,Rv)
    print(f"{name:10s} {k:14s} left={s['left']:.2f} mirror-matches={nm} on-plane={npl} oddfrac={odd_fraction(m,p):.3f} closest={np.round(best[2:],3)}")

import json, numpy as np, time
from solve import *
from lcdm import sigma_G
out={}
for H in [0.3,0.4]:
    g=Grid(H)
    for s in ([0,1,2,4] if H==0.3 else [2]):
        t=time.time(); S=Solver(g,s); curve=[]
        for lg in np.arange(-14,8.01,0.25):
            sc=S.score(S.field(10**lg)); sc["lam"]=10**lg; curve.append(sc)
        d,best=S.to_target(1.0)
        key=f"H{H}_s{s}"; out[key]=dict(H=H,s=s,curve=curve,at1x=best,sigma_LCDM_gauss=float(sigma_G(s)) if s>0 else float(sigma_G(H)))
        if H==0.3: np.save(f"delta_{key}.npy",d)
        print(key,{k:round(v,4) if isinstance(v,float) else v for k,v in best.items()},"minleft",round(min(c["left"] for c in curve),4),f"{time.time()-t:.0f}s",flush=True)
json.dump(out,open("tradeoff.json","w"))

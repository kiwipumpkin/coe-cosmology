import json, time, numpy as np
from field_multi import *
S=2.0; res={}
def mk(where,half):
    gs=[]
    if where in("inside","both"): gs.append(G(0.3,0.5,13.8,half))
    if where in("outside","both"): gs.append(G(0.6,14.8,30.0,half))
    return gs
for where in ["inside","outside","both"]:
    for half in [0,1,-1]:
        t=time.time(); gs=mk(where,half); r=solve(gs,S)
        key=f"{where}_{'mirror' if half==0 else ('onesided+' if half>0 else 'onesided-')}"
        if r is None: res[key]=dict(reach=False); print(key,"cannot reach 1x",flush=True); continue
        ds,L=r; st=stats(gs,ds); res[key]=dict(reach=True,left=L,parts=st,M_plus=sum(p["M_plus"] for p in st),M_minus=sum(p["M_minus"] for p in st))
        np.savez(f"field_{key}.npz",*ds,*[g.P for g in gs])
        print(key,f"left={L:.3f} M+={res[key]['M_plus']:.2e} M-={res[key]['M_minus']:.2e}",[{k:(round(v,4) if abs(v)<10 else f'{v:.2e}') for k,v in p.items()} for p in st],f"{time.time()-t:.0f}s",flush=True)
json.dump(res,open("fields.json","w"),indent=1)

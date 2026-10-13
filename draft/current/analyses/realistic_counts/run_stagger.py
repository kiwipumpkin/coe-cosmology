import json, numpy as np
from field_multi import *
res=json.load(open("fields.json"))
for name,hin,hout in [("both_staggered_in+out-",1,-1),("both_staggered_in-out+",-1,1)]:
    gs=[G(0.3,0.5,13.8,hin),G(0.6,14.8,30.0,hout)]; r=solve(gs,2.0)
    if r is None: res[name]=dict(reach=False); print(name,"cannot"); continue
    ds,L=r; st=stats(gs,ds); res[name]=dict(reach=True,left=L,parts=st,M_plus=sum(p["M_plus"] for p in st),M_minus=sum(p["M_minus"] for p in st))
    np.savez(f"field_{name}.npz",*ds,*[g.P for g in gs]); print(name,res[name]["M_plus"],[ (p["rms"],p["peak"]) for p in st],flush=True)
json.dump(res,open("fields.json","w"),indent=1)

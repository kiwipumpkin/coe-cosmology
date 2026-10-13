import json, numpy as np
from field_multi import *
res=json.load(open("fields.json"))
for name,a,b in [("inside_staggered_near+far-",1,-1),("inside_staggered_near-far+",-1,1)]:
    gs=[G(0.3,11.5,13.8,a),G(0.3,0.5,9.8,b)]   # one side only near the ring (>=11.5 Gpc), the other side only <=9.8 Gpc (>15% closer)
    r=solve(gs,2.0)
    if r is None: res[name]=dict(reach=False); print(name,"cannot"); continue
    ds,L=r; st=stats(gs,ds); res[name]=dict(reach=True,left=L,parts=st,M_plus=sum(p["M_plus"] for p in st),M_minus=sum(p["M_minus"] for p in st))
    np.savez(f"field_{name}.npz",*ds,*[g.P for g in gs]); print(name,f"{res[name]['M_plus']:.2e}",[(round(p["rms"],4),round(p["peak"],3)) for p in st],flush=True)
json.dump(res,open("fields.json","w"),indent=1)

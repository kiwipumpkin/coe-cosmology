import json,numpy as np
from cb_core import *
res=[]
for mdl in ["chess","einasto"]:
    for DT,h,nth in [(0.5,LAM/6,360),(1.0,LAM/6,360),(1.5,LAM/6,360),(1.0,LAM/9,360),(1.0,LAM/6,3600)]:
        w,_=pattern(mdl,0.10,DT=DT,h=h,nth=nth)
        r=dict(model=mdl,DT=DT,h_Mpc=round(h*1000,1),nth=nth,own=rms(w)/FLOOR,m1_12=band(w,1,12),m13_180=band(w,13,180),m13up=band(w,13,99999),sym=sym_share(w))
        print(r,flush=True); res.append(r)
json.dump(res,open("convergence.json","w"),indent=1)

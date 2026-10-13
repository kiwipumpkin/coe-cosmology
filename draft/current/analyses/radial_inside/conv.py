import json
from ri_core import *
out=[]
for lam in [LAM,0.150]:
  for DT,h,nth in [(1.0,LAM/6,360),(2.0,LAM/6,360),(1.0,LAM/9,360),(1.0,LAM/6,3600)]:
    w,wus=pattern_r(0.10,lam,DT=DT,h=h,nth=nth)
    r=dict(lam=lam,DT=DT,h=h,nth=nth,own=rms(w)/FLOOR,m1_12=band(w,1,12),m1_5=band(w,1,5),m13_180=band(w,13,180),m13up=band(w,13,99999),sym=sym_share(w),sym_1_12=sym_share(w,1,12),us=rms(wus-wus.mean())/FLOOR)
    print(r,flush=True); out.append(r)
json.dump(out,open("convergence.json","w"),indent=1,default=float)

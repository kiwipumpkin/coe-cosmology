import time
from ps_core import *
gs=CONFIGS["in10"]()
for B in [0.02,0.03,0.035,0.05]:
    t=time.time(); ds,lp,ld=solve_dual(gs,B); d=ds[0]
    print(B,round(lp,3),round(ld,3),"peak",abs(d).max(),"frac at bound",(abs(d)>0.98*B).mean(),round(time.time()-t,1),"s",flush=True)

import json,time,numpy as np,half_fit as H
from exact2 import node
p=json.load(open('half_ckpt_r2free.json'))['best']['p']
for L in (1500,2500):
  for dk in (0.5,0.25,0.12):t=time.time();print(L,dk,round(node((p,L,dk)),4),round(time.time()-t,1),flush=True)

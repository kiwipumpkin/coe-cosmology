import time; from bound_core import *
t=time.time()
for mc in (0.02,0.3,1,2,5):
  for DR in (1.5,2.7,5,9):
    for LD in (0.07,0.3,1,3):
      rr=run(DR,LD,mc); print(mc,DR,LD,"vrmin %.3f vesc2 %.3f"%(rr["vrmin"],rr["vesc2"]), "R2 %.3f"%score(dict(DR=DR,LD=LD,mc=mc),"mig")[0], "R2n %.3f"%score(dict(DR=DR,LD=LD,mc=mc),"mig","mean")[0])
print(time.time()-t)

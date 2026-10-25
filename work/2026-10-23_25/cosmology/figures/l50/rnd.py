import numpy as np
src=open('jk.py').read().split("allk=np.ones")[0];exec(src)
allk=np.ones(len(z),bool);S=[np.where(grid==c)[0][0] for c in (295,300,305,310,315)]
idx=np.where(np.abs((l-300+180)%360-180)<=65)[0];rng=np.random.default_rng(1);v=[]
for t in range(100):
  k=allk.copy();k[rng.choice(idx,10,replace=False)]=False;A1,B1=curves(k);v.append((B1-BLm)[S].mean())
print('random10 B %.4f +- %.4f'%(np.mean(v),np.std(v)))

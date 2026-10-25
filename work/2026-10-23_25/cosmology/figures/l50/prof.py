import numpy as np
src=open('jk.py').read().split("allk=np.ones")[0];exec(src)
allk=np.ones(len(z),bool);A0,B0=curves(allk);rA=A0-AFm;rB=B0-BLm
for g,a,bb in zip(grid,rA,rB):
  if 240<=g<=350: print(g,'%+.4f %+.4f'%(a,bb))

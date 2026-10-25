import time,numpy as np; from lc import *
z=[-1,-1,0,0.3,1.0,0.5]  # exp law guess
for zz in ([-1,-1,0,np.log10(1.5),1.0,-0.3],[-0.5,-0.5,0.3,np.log10(2),1.0,0.0]):
  L,tn=make('exp',zz); o=obs(L,tn); t0=time.time(); af=migkas(L,tn); print('U',o['U'],'H',o['Hpar'],o['Hper'], 'time',time.time()-t0, None if af is None else F.chi2_l0(af))

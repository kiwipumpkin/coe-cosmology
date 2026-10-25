import time,numpy as np; from lc import *
L,tn=make('exp',[-1,-1,0,np.log10(1.5),1.0,-0.3])
for th in (0,90,180):
  for tm in (1e-2,1e-4):
    t0=time.time(); r=ray(L,tn,np.radians(th),tmin=tm*tn,rtol=1e-8); lr=np.log(L.rho_rel(r['t'],r['x'])/obs(L,tn)['rho'])
    print(th,tm,'%.1fs'%(time.time()-t0),r['stop'],len(r['t']),r['t'][-1]/tn,int(r['nc'].max()),'lnz1 end',r['lnz1'][-1],'lnrho end',lr[-1],flush=True)

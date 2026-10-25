import time,json,numpy as np; from lc import *
L,tn=make('exp',[-1,-1,0,np.log10(1.5),1.0,-0.3]); t0=time.time()
res,RD=cmb(L,tn,np.arange(0,360,30.),tmode='A'); print('t',time.time()-t0,'K',res.get('K'),res.get('tstar'),tn)
for r in RD: f=np.isfinite(r['lJ']); print(round(np.degrees(r['th'])),r['stop'],'max lJ',r['lJ'][f].max() if f.any() else None)
if res.get('K'):
  for q in res['rows']: print({k:v for k,v in q.items()})
  print(score(res['rows'],180,None))
else: print(res.get('K_scan','')[:5])
import pickle; pickle.dump(RD,open('t_rd.pkl','wb'))

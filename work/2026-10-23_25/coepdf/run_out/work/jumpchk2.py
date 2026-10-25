import json,time,numpy as np,half_fit as H,acoustic_model as M
C=json.load(open('half_ckpt_r2free.json'));p=C['best']['p']
t=M.EX.terms(p);Dex=t['scalar']+t['mlens'];L=np.array([30,40,50])
print('Fwd exact',Dex[L-2] if len(Dex)==49 else Dex[L], len(Dex))
print('limber*ratio',M.limber(p,L)*np.array([0.6284,0.7012,0.7517]))
t0=time.time();n,r=M.exact_nodes(p,nodes=(1000,1500,2000,2500));print(time.time()-t0,dict(zip(n,np.round(r,4))))

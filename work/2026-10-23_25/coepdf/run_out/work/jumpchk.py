import json,time,numpy as np,half_fit as H,acoustic_model as M
C=json.load(open('half_ckpt_r2free.json'));p=C['best']['p'];H.CORR[:]=C['corr']
D=M.Dfield(p)*H.CORR+M.Dfix;print('stale jump 50->51',D[49:53])
t=time.time();n,r=M.exact_nodes(p,nodes=(30,40,50,51,60,75,100,150,200,300,400,500,700));print('t',time.time()-t);print(dict(zip(n,np.round(r,4))))

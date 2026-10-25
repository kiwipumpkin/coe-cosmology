import json,numpy as np;from jac import *
P=json.load(open('inv1_plan.json'));out={}
for U in P:
    for n,v in P[U].items():
        if not n.startswith(U+'_inv'):continue
        f=f'runs/{n}.npz'
        if not os.path.exists(f):out[n]='run missing/failed';continue
        am,ag=model(f,v['C'],v['DRV']);i1,c1,_=bestaxis(am,ag,v['Ag']);i2,c2,A2=bestaxis(am,ag)
        out[n]=dict(params=v,chi2_planAg=c1,axis_planAg=[float(x) for x in hp.pix2ang(NA,i1,lonlat=True)],chi2_Agfree=c2,Ag_free=A2,axis_Agfree=[float(x) for x in hp.pix2ang(NA,i2,lonlat=True)],a_motion=am.tolist(),a_grav=ag.tolist(),Dl_motion=Dl_al0(am).tolist())
        print(n,c1,c2,A2,flush=True)
json.dump(out,open('inv1_results.json','w'),indent=1)

import json,numpy as np;from scipy.optimize import lsq_linear
from jac import *
R=json.load(open('jac_results.json'));out={};jobs=[]
for U,b in BASE.items():
    J=np.array(R[U]['J_muK']);Jw=J/sig[:,None]
    am0,ag0=model(b['run'],b['C'],b['DRV']);ia,c0,_=bestaxis(am0,ag0,1.0);r=(Ga[L2,ia]-(am0+ag0)[L2])/sig
    lf0=np.log10(b['frac']);lo=[-.05]*4+[max(1-lf0,-1),1.5-b['DRV'],-1e4];hi=[.05]*4+[min(5-lf0,1),6-b['DRV'],1e4]
    s=lsq_linear(Jw,r,bounds=(lo,hi));d=s.x
    out[U]=dict(bounds=[lo,hi],step=d.tolist(),lin_pred_chi2=float(((r-Jw@d)**2).sum()),start_chi2=c0,active=s.active_mask.tolist())
    for t in(1.0,0.5,0.25):
        C=[b['C'][i]+t*d[i] for i in range(4)];fr=10**(lf0+t*d[4]);dr=b['DRV']+t*d[5]
        n=f'{U}_inv_t{t}';jobs.append(f'{n} {C[0]:.6g} {C[1]:.6g} {C[2]:.6g} {C[3]:.6g} {fr:.6g} {dr:.6g}');out[U][n]=dict(C=C,frac=fr,DRV=dr,Ag=1+t*d[6])
    print(U,out[U])
json.dump(out,open('inv1_plan.json','w'),indent=1);open('jobs_inv1.txt','w').write('\n'.join(jobs)+'\n')

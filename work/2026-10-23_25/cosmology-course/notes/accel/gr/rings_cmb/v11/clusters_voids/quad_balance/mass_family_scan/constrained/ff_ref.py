import json,numpy as np
from ff_prof import fit
from multiprocessing import Pool
if __name__=='__main__':
    J=[(r,-10,100) for r in np.linspace(1.045,1.075,13)]+[(r,0,360) for r in np.r_[np.linspace(0.6,1.0,9),np.linspace(1.045,1.085,9)]]
    with Pool(4) as P:out=P.map(fit,J)
    R=json.load(open('fullfree_results.json'));R['refine']=out;json.dump(R,open('fullfree_results.json','w'),indent=1)
    for d in out:print(d['l0range'][1],'%.4f'%d['ratio'],'%.5f'%d['chi2'],'%.4f'%d['R2'],np.round(d['theta'],3).tolist())

import os,sys,json,numpy as np,warnings;warnings.filterwarnings('ignore')
H=os.getcwd();os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/density')
exec(open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/migkas_prof.py').read().split('rng=np.random')[0]);os.chdir(H)
from multiprocessing import Pool
def fit(a):
    rat,l0lo,l0hi=a;lo=LO.copy();hi=HI.copy();lo[1]=l0lo;hi[1]=l0hi
    f=lambda th:1e9 if np.any(th<lo) or np.any(th>hi) else np.sum(((model(th,rat)-Y)/S)**2)
    rng=np.random.default_rng(1);best=None
    for k in range(40):
        r=minimize(f,lo+(hi-lo)*rng.random(6),method='Nelder-Mead',options=dict(maxiter=6000,xatol=1e-8,fatol=1e-8))
        if best is None or r.fun<best.fun:best=r
    p=best.x
    for j in range(3):p=minimize(f,p,method='Nelder-Mead',options=dict(maxiter=6000,xatol=1e-10,fatol=1e-10)).x
    m=model(p,rat);return dict(ratio=rat,l0range=[l0lo,l0hi],chi2=float(f(p)),R2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)),theta=p.tolist())
if __name__=='__main__':
    rats=sorted(set(np.r_[1+np.geomspace(1e-6,0.3,30),np.linspace(1.3,2.2,10),[1.0,0.99,0.95,0.9]].round(9)))
    jobs=[(r,a,b) for r in rats for a,b in [(-10,100),(0,360)]]
    with Pool(4) as P:out=P.map(fit,jobs)
    json.dump(dict(profile=out),open('fullfree_results.json','w'),indent=1)
    for d in out:print(d['l0range'][1],'%.7f'%d['ratio'],'%.4f'%d['chi2'],'%.4f'%d['R2'],np.round(d['theta'],3).tolist())

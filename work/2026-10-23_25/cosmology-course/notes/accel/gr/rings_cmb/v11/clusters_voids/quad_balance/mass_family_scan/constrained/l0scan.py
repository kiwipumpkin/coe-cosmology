import os,json,numpy as np,warnings;warnings.filterwarnings('ignore')
H=os.getcwd();os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/density')
exec(open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/migkas_prof.py').read().split('rng=np.random')[0]);os.chdir(H)
from multiprocessing import Pool
R0=json.load(open('interp_results.json'));L0S=json.loads(os.environ['L0S']) if 'L0S' in os.environ else list(range(15,101,5))
idx=[0,2,3,4,5]
def fit(rat,l0,n=40):
    rng=np.random.default_rng(12345);f=lambda p:chi(np.insert(p,1,l0),rat);best=None
    for k in range(n):
        th=LO+(HI-LO)*rng.random(6);r=minimize(f,th[idx],method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-7,fatol=1e-7))
        if best is None or r.fun<best.fun:best=r
    p=best.x
    for j in range(3):p=minimize(f,p,method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-9,fatol=1e-9)).x
    th=np.insert(p,1,l0);m=model(th,rat)
    return float(chi(th,rat)),float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)),th.tolist()
def work(a):
    k,l0=a;c,r2,th=fit(R0[k]['ratio'],l0);return dict(run=k,l0=l0,ratio=R0[k]['ratio'],chi2=c,R2=r2,theta=th,theta_names=['log10u','l0','bf','vp','lpp','sb'])
if __name__=='__main__':
    out='l0scan_results.json';D=json.load(open(out)) if os.path.exists(out) else []
    done={(d['run'],d['l0']) for d in D}
    jobs=[(k,l) for k in R0 if k!='_meta' and 'ratio' in R0[k] for l in L0S if (k,l) not in done]
    with Pool(7) as p:
        for i,d in enumerate(p.imap_unordered(work,jobs)):
            D.append(d)
            if i%20==0 or i==len(jobs)-1:json.dump(D,open(out,'w'))
    json.dump(D,open(out,'w'));print('done',len(D))

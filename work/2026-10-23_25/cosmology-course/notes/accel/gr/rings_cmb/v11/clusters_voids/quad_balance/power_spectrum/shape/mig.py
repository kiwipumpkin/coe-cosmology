# Migkas chi2/R2 for given H_par/H_perp ratios; identical model/data/fitter as mass_family_scan/migkas_prof.py
import os,sys,json,numpy as np,warnings;warnings.filterwarnings('ignore')
rats=[float(a) for a in sys.argv[2:]];outf=os.path.abspath(sys.argv[1])
os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/density')
src=open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/migkas_prof.py').read()
exec(src.split('rng=np.random')[0])
rng=np.random.default_rng(5);out={}
for rat in rats:
    best=None
    for k in range(25):
        r=minimize(chi,LO+(HI-LO)*rng.random(6),args=(rat,),method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-7,fatol=1e-7))
        if best is None or r.fun<best.fun:best=r
    m=model(best.x,rat);out[rat]=dict(chi2=float(best.fun),R2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)));print(rat,out[rat],flush=True)
json.dump(out,open(outf,'w'),indent=1)

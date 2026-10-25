# Stable Migkas profile chi2(ratio): identical model/data/bounds/fitter as mass_family_scan/migkas_prof.py; more restarts + fixed per-ratio seed + Nelder-Mead polish.
import os,sys,json,numpy as np,warnings;warnings.filterwarnings('ignore')
outf=os.path.abspath(sys.argv[1]);rats=[float(a) for a in sys.argv[2:]];NR=int(os.environ.get('NR','150'))
os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/density')
src=open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/migkas_prof.py').read()
exec(src.split('rng=np.random')[0]);out={}
for rat in rats:
    rng=np.random.default_rng(12345);fs=[]
    for k in range(NR):
        r=minimize(chi,LO+(HI-LO)*rng.random(6),args=(rat,),method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-7,fatol=1e-7))
        for j in range(3): r=minimize(chi,r.x,args=(rat,),method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-8,fatol=1e-8))
        fs.append(r)
    fs.sort(key=lambda r:r.fun);b=fs[0];m=model(b.x,rat)
    out[rat]=dict(chi2=float(b.fun),R2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)),l0=float(b.x[1]),theta=b.x.tolist(),n_within_0p1=int(sum(f.fun<b.fun+0.1 for f in fs)),N=len(Y))
    print(rat,out[rat]['chi2'],out[rat]['l0'],out[rat]['n_within_0p1'],flush=True)
json.dump(out,open(outf,'w'),indent=1)

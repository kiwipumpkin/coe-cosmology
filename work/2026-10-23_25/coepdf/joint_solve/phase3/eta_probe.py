# NEW: profile chi2 vs eta below the 0.03 guard of cw.sky_locked (guard relaxed in-memory only; source file untouched). Squared model A=af^2/<af^2>, L locked.
import numpy as np,json,types
from scipy.optimize import minimize
src=open('/workspace/cosmology-course/notes/catwise-eta-density-migkas.py').read().replace('0.03 < eta < 8.0','1e-4 < eta < 8.0')
cw=types.ModuleType('cw');cw.__file__='/workspace/cosmology-course/notes/catwise-eta-density-migkas.py';exec(compile(src,cw.__file__,'exec'),cw.__dict__)
B=json.load(open('/workspace/coepdf/joint_solve/phase3/migkas_refit_sq.json'))['refit_sq'];DL=cw.RV
def chi2(x,eta):
    D,mc,l0=10**x[0],10**x[1],x[2];o=cw.sky_locked(D,eta,mc,DL)
    if o is None:return 800.
    af=o[0];a2=af**2/np.mean(af**2);return float(np.sum(cw.W*(cw.shifted(a2,l0)[cw.IMIG]-cw.Y)**2))
x=[np.log10(B['D']),np.log10(B['mc']),B['l0']];out=[]
for eta in [0.03001,0.02,0.01,0.005,0.002,0.001]:
    r=min((minimize(chi2,x0,args=(eta,),method='Nelder-Mead',options=dict(maxfev=500,xatol=1e-4,fatol=1e-4)) for x0 in (x,[x[0],x[1]+.3,x[2]],[x[0]+.1,x[1],x[2]])),key=lambda r:r.fun)
    x=list(r.x);out.append(dict(eta=eta,chi2=r.fun,D=10**x[0],mc=10**x[1],l0=x[2]));print(out[-1],flush=True)
json.dump(out,open('/workspace/coepdf/joint_solve/phase3/eta_probe.json','w'),indent=1)

# Migkas chi2 profile vs the H_par/H_perp ratio at our position (the only way mass enters the Migkas model), other params free.
import os,json,numpy as np,warnings;warnings.filterwarnings('ignore')
os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/density')
exec(open('../phys/phys.py').read().split('def age')[0])
from scipy.optimize import minimize
C=299792.458;LO=np.array([np.log10(0.3),-60,-0.5,0,0,-1]);HI=np.array([np.log10(30),100,0.5,3000,360,1])
def model(th,rat):
    lu,l0,bf,vp,lpp,sb=th;u=10**lu;G=rest_grad(rat*u,u,bf,1.0);Hpa,Hpe=G[0,0],G[1,1]
    t=np.radians(L-l0);tt=np.radians(np.arange(360)-l0);Hb=np.mean(Hpe+(Hpa-Hpe)*np.cos(tt)**2)
    return ((Hpe+(Hpa-Hpe)*np.cos(t)**2)/Hb*(1+bf*WI*np.cos(t))-vp/C*WI*np.sqrt(1-sb**2)*np.cos(np.radians(L-lpp)))**2
def chi(th,rat):
    if np.any(th<LO) or np.any(th>HI):return 1e9
    return np.sum(((model(th,rat)-Y)/S)**2)
rng=np.random.default_rng(5);out={}
for rat in [1.0,1.005,1.01,1.02,1.03,1.04,1.05,1.06,1.07,1.08,1.1,1.13,1.16,1.2]:
    best=None
    for k in range(25):
        r=minimize(chi,LO+(HI-LO)*rng.random(6),args=(rat,),method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-7,fatol=1e-7))
        if best is None or r.fun<best.fun:best=r
    m=model(best.x,rat);out[rat]=dict(chi2=float(best.fun),R2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)));print(rat,out[rat],flush=True)
json.dump(out,open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/migkas_vs_ratio.json','w'),indent=1)

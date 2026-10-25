import json,numpy as np,warnings,emcee;warnings.filterwarnings('ignore')
exec(open('../phys/phys.py').read().split('fs=np.linspace')[0])
OS=0.0
FG=np.linspace(0,1.5,151);RG=np.array([ratio(f,1.0) for f in FG]);np.save('ratio_f.npy',np.c_[FG,RG])
C=299792.458
def model(th):
    f,lu,l0,bf,vp,lp,sb=th; u=10**lu; rat=np.interp(f,FG,RG)
    G=rest_grad(rat*u,u,bf,1.0);Hpa,Hpe=G[0,0],G[1,1]
    t=np.radians(L-l0);tt=np.radians(np.arange(360)-l0)
    Hb=np.mean(Hpe+(Hpa-Hpe)*np.cos(tt)**2);H=Hpe+(Hpa-Hpe)*np.cos(t)**2
    cb=np.sqrt(1-sb**2)
    return (H/Hb*(1+bf*WI*np.cos(t))-vp/C*WI*cb*np.cos(np.radians(L-lp)))**2
LO=np.array([0,np.log10(0.3),-60,-0.5,0,0,-1]);HI=np.array([1.5,np.log10(30),100,0.5,3000,360,1])
def lnp(th):
    if np.any(th<LO) or np.any(th>HI): return -np.inf
    return -0.5*np.sum(((model(th)-Y)/S)**2)
# multistart optimiser
from scipy.optimize import minimize
rng=np.random.default_rng(1);best=None
for k in range(60):
    x0=LO+(HI-LO)*rng.random(7)
    r=minimize(lambda x:-lnp(x) if np.isfinite(lnp(x)) else 1e9,x0,method='Nelder-Mead',options=dict(maxiter=6000,xatol=1e-6,fatol=1e-6))
    if best is None or r.fun<best.fun: best=r
xb=best.x;print('opt',xb,2*best.fun)
nw=40;p0=xb+1e-3*(HI-LO)*rng.standard_normal((nw,7));p0=np.clip(p0,LO+1e-6,HI-1e-6)
sm=emcee.EnsembleSampler(nw,7,lnp);sm.run_mcmc(p0,12000,progress=False)
ch=sm.get_chain(discard=3000,thin=10,flat=True);np.save('chain.npy',ch)
try: tau=sm.get_autocorr_time(discard=3000,quiet=True).tolist()
except Exception as e: tau=str(e)
names=['core_frac','log10_uD','l0','beta_f','v_pec_kms','l_pec','sin_b_pec']
q=lambda a:np.percentile(a,[16,50,84]).tolist()
out=dict(best=dict(zip(names,xb)),chi2=2*best.fun,post={n:q(ch[:,i]) for i,n in enumerate(names)},tau=tau)
m=model(xb);out['R2']=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2))
vr=ch[:,4]*np.sqrt(1-ch[:,6]**2);out['post']['v_pec_ring_kms']=q(vr)
out['post']['ratio']=q(np.interp(ch[:,0],FG,RG))
out['corr']=np.corrcoef(ch.T).round(2).tolist()
json.dump(out,open('point_fit.json','w'),indent=1,default=float);print(json.dumps(out,indent=1,default=float))

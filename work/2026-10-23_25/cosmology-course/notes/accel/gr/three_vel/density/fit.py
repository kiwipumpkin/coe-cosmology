import json,sys,pickle,numpy as np,warnings,emcee;warnings.filterwarnings('ignore')
from scipy.interpolate import RegularGridInterpolator as RGI,interp1d
exec(open('../phys/phys.py').read().split('def age')[0])
g=pickle.load(open('grid.pkl','rb'));DG,PG,RR=g['DG'],g['PG'],g['R'];C=299792.458
def rat_fn(k):
    if k=='exp':
        f=interp1d(np.log(DG),RR[k][:,0,0]);return lambda D,p:float(f(np.log(D)))
    Z=RR[k][...,0];f=RGI((np.log(DG),np.log(PG[k])),Z,bounds_error=False,fill_value=np.nan)
    return lambda D,p:float(f([[np.log(D),np.log(p)]])[0])
def make(k):
    rf=rat_fn(k);pl,ph=(0,0) if k=='exp' else (np.log(PG[k][0]),np.log(PG[k][-1]))
    LO=np.array([np.log(15),pl,np.log10(0.3),-60,-0.5,0,0,-1]);HI=np.array([np.log(400),ph,np.log10(30),100,0.5,3000,360,1])
    def model(th):
        lD,lp_,lu,l0,bf,vp,lpp,sb=th;rat=rf(np.exp(lD),np.exp(lp_))
        if not np.isfinite(rat):return None
        u=10**lu;G=rest_grad(rat*u,u,bf,1.0);Hpa,Hpe=G[0,0],G[1,1]
        t=np.radians(L-l0);tt=np.radians(np.arange(360)-l0);Hb=np.mean(Hpe+(Hpa-Hpe)*np.cos(tt)**2)
        return ((Hpe+(Hpa-Hpe)*np.cos(t)**2)/Hb*(1+bf*WI*np.cos(t))-vp/C*WI*np.sqrt(1-sb**2)*np.cos(np.radians(L-lpp)))**2
    def lnp(th):
        if k=='exp':th=th.copy();th[1]=0
        if np.any(th<LO) or np.any(th>HI):return -np.inf
        m=model(th);return -np.inf if m is None else -0.5*np.sum(((m-Y)/S)**2)
    return model,lnp,LO,HI
if __name__=='__main__':
    from scipy.optimize import minimize
    k=sys.argv[1];model,lnp,LO,HI=make(k);rng=np.random.default_rng(3);best=None
    for i in range(60):
        x0=LO+(HI-LO)*rng.random(8)
        r=minimize(lambda x:-lnp(x) if np.isfinite(lnp(x)) else 1e9,x0,method='Nelder-Mead',options=dict(maxiter=8000,xatol=1e-7,fatol=1e-7))
        if best is None or r.fun<best.fun:best=r
    xb=best.x;nw=40;p0=np.clip(xb+1e-3*(HI-LO+1e-9)*rng.standard_normal((nw,8)),LO+1e-7,HI-1e-7)
    sm=emcee.EnsembleSampler(nw,8,lnp);sm.run_mcmc(p0,8000,progress=False);ch=sm.get_chain(discard=2000,thin=10,flat=True);np.save(f'chain_{k}.npy',ch)
    m=model(xb);q=lambda a:np.percentile(a,[16,50,84]).tolist()
    nm=['ln_D_Gpc','ln_shape','log10_uD','l0','beta_f','v_pec_kms','l_pec','sin_b_pec']
    rat=rat_fn(k)(np.exp(xb[0]),np.exp(xb[1]))
    out=dict(law=k,best=dict(zip(nm,xb)),D_Gpc=float(np.exp(xb[0])),shape=float(np.exp(xb[1])),ratio=rat,chi2=2*best.fun,R2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)),
      post={n:q(ch[:,i]) for i,n in enumerate(nm)},post_D=q(np.exp(ch[:,0])),post_shape=q(np.exp(ch[:,1])),vpec_at_cap=bool(xb[5]>2990),frac_chain_vp_gt2900=float(np.mean(ch[:,5]>2900)))
    json.dump(out,open(f'fit_{k}.json','w'),indent=1,default=float);print(k,out['chi2'],out['R2'],out['D_Gpc'],out['shape'],rat,xb[3:],out['vpec_at_cap'])

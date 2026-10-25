import json,numpy as np,warnings,emcee;warnings.filterwarnings('ignore')
exec(open('../phys/phys.py').read().split('fs=np.linspace')[0])
OS=0.0
FG=np.linspace(0,1.5,151);RG=np.array([ratio(f,1.0) for f in FG]);
C=299792.458
def model(th):
    f,lu,l0,bf,vp,lp,sb=th; u=10**lu; rat=np.interp(f,FG,RG)
    G=rest_grad(rat*u,u,bf,1.0);Hpa,Hpe=G[0,0],G[1,1]
    t=np.radians(L-l0);tt=np.radians(np.arange(360)-l0)
    Hb=np.mean(Hpe+(Hpa-Hpe)*np.cos(tt)**2);H=Hpe+(Hpa-Hpe)*np.cos(t)**2
    cb=np.sqrt(1-sb**2)
    return (H/Hb*(1+bf*WI*np.cos(t))-vp/C*WI*cb*np.cos(np.radians(L-lp)))**2
LO=np.array([0,np.log10(4.0),-60,-0.5,0,0,-1]);HI=np.array([1.5,np.log10(30),100,0.5,3000,360,1])
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

from scipy.stats import chi2 as C2,norm
xu=np.array([json.load(open('point_fit.json'))['best'][n] for n in ['core_frac','log10_uD','l0','beta_f','v_pec_kms','l_pec','sin_b_pec']])
wm=np.average(Y,weights=1/S**2);cf=float(np.sum(((Y-wm)/S)**2))
def st(x):
    m=model(x);c=float(np.sum(((m-Y)/S)**2));R2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2));d=cf-c;sig=float(norm.isf(C2.sf(d,6)/2)) if d>0 else 0.0
    return dict(chi2=c,R2=R2,dchi2_vs_flat=d,sigma_vs_flat=sig)
out=dict(names=['core_frac','log10_uD','l0','beta_f','v_pec_kms','l_pec','sin_b_pec'],constrained=dict(x=xb.tolist(),u=10**xb[1],**st(xb)),unconstrained=dict(x=xu.tolist(),u=10**xu[1],**st(xu)),flat_chi2=cf,n=len(Y),
  sigma_note='flat = weighted-mean constant; delta chi2 vs flat with 6 extra params (7 minus the constant), converted to two-sided Gaussian sigma')
json.dump(out,open('point_fit_c4.json','w'),indent=1,default=float);print(json.dumps(out,indent=1,default=float))
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
o=np.argsort(L);fig,ax=plt.subplots(figsize=(13,8))
ax.errorbar(L,Y,yerr=S,fmt='o',color='k',ms=4,alpha=.6,label=f'Migkas A/A_all data ({len(Y)} bins)')
lg=np.linspace(L.min(),L.max(),721)
for k,col,lab in (('unconstrained','C0','unconstrained best fit'),('constrained','C3','constrained: our speed at least 4 c')):
    d=out[k];x=np.array(d['x'])
    ax.plot(L[o],model(x)[o],'-',color=col,lw=2.5,label=f"{lab}: speed {d['u']:.3f} c, core mass {x[0]:.4f} of critical, R² {d['R2']:.3f}, χ² {d['chi2']:.2f}, {d['sigma_vs_flat']:.1f}σ vs flat")
ax.grid(alpha=.3);ax.set_xlabel('Galactic longitude of cluster bin (deg)');ax.set_ylabel('A / A_all')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10,frameon=False)
fig.savefig('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/migkas_c4_fit.png',dpi=110,bbox_inches='tight')

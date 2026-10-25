import json,numpy as np,warnings;warnings.filterwarnings('ignore')
src=open('point.py').read().split('# multistart')[0].replace("FG=np.linspace(0,1.5,151);RG=np.array([ratio(f,1.0) for f in FG]);np.save('ratio_f.npy',np.c_[FG,RG])","FG,RG=np.load('ratio_f.npy').T")
exec(src)
from scipy.optimize import minimize
d=json.load(open('point_fit.json'));xb=np.array(list(d['best'].values()))
def fitfix(fix):
    free=[i for i in range(7) if i not in fix];best=None
    rng=np.random.default_rng(2)
    for k in range(25):
        x0=xb.copy()
        for i in free: x0[i]=LO[i]+(HI[i]-LO[i])*rng.random() if k else xb[i]
        for i,v in fix.items(): x0[i]=v
        def fn(z):
            x=x0.copy();x[free]=z;v=lnp(x);return -v if np.isfinite(v) else 1e9
        r=minimize(fn,x0[free],method='Nelder-Mead',options=dict(maxiter=6000))
        if best is None or r.fun<best.fun: best=r;bx=x0.copy();bx[free]=r.x
    return bx
st=lambda m:(float(np.sum(((m-Y)/S)**2)),float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)))
cases={'best fit, all free':xb,'own motion fixed 370 km/s to Planck apex':fitfix({4:370,5:264.021,6:np.sin(np.radians(48.253))}),
 'frame speed 0, own motion free':fitfix({3:0.0})}
res={k:dict(chi2_R2=st(model(v)),params=v.tolist()) for k,v in cases.items()};json.dump(res,open('cases.json','w'),indent=1);print(json.dumps(res,indent=1))
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':16})
fig,ax=plt.subplots(figsize=(13,8));o=np.argsort(L);ax.errorbar(L,Y,yerr=[MLO,MHI],fmt='o',color='k',ms=4,label='Migkas A/A_all (71 bins)')
for (k,v),c in zip(cases.items(),['C0','C2','C3']):
    a,b=st(model(v));ax.plot(L[o],model(v)[o],c,lw=3,label=f'{k}: χ²={a:.1f}, R²={b:.2f}')
ax.set_xlabel('Galactic longitude (deg)');ax.set_ylabel('A / A_all');ax.legend(loc='lower center',bbox_to_anchor=(.5,1),frameon=False);fig.tight_layout();fig.savefig('point_fit_migkas.png',dpi=120)
import corner;ch=np.load('chain.npy')
lab=['core mass /\ncritical','log10 space speed\nat D (c)','axis\nlongitude','frame speed\n(c)','own speed\n(km/s)','own-motion\nlongitude','sin(own-motion\nlatitude)']
fg=corner.corner(ch,labels=lab,quantiles=[.16,.5,.84],show_titles=True,title_fmt='.3g',label_kwargs=dict(fontsize=12));fg.savefig('point_corner.png',dpi=90)

import sys,json,numpy as np,healpy as hp;sys.path.insert(0,'../../rotfit')
import fit
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
o=np.load('../b_1100_k0.npz')['o'];fit.T1=CubicSpline(o[:,0],1/o[:,1])
vec,lb,ang,dip,sky,T0=fit.vec,fit.lb,fit.ang,fit.dip,fit.sky,fit.T0
aC=vec(0,15);CoE=-aC;PH=vec(264.02,48.25)
def rep(v,name):
    m=sky(aC,v);d=dip(m);return dict(case=name,motion_lb=None if v is None else lb(v),amp_xPlanck=float(np.linalg.norm(d)*T0/3362),hot_lb=lb(d),offset_from_Planck_hot_deg=float(ang(d,PH)),offset_from_CoE_deg=float(ang(d,CoE)))
dm=dip(sky(aC,None))
# motion-only dipole: isotropic sky + motion
r=[rep(None,'model only')]
r.append(dict(case='motion only (toward Planck hot)',amp_xPlanck=float(fit.beta*T0/3362),hot_lb=lb(PH)))
r+= [rep(CoE,'motion toward CoE (aligned)'),rep(-CoE,'motion away from CoE (opposite)'),rep(PH,'motion toward observed Planck dipole (264.02, +48.25)')]
b=min((minimize(lambda x:np.linalg.norm(dip(sky(aC,vec(*x)))),s,method='Nelder-Mead',options=dict(xatol=.01,fatol=1e-9)) for s in [[0,15],[10,10],[350,20]]),key=lambda q:q.fun)
r.append(rep(vec(*b.x),'motion chosen to minimise (2-angle fit)'))
for x in r:print(x)
json.dump(r,open('motion_results.json','w'),indent=1,default=float)
# plot dT/T along great circle through CoE and Planck hot
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
e1=CoE;e2=PH-(PH@e1)*e1;e2/=np.linalg.norm(e2);c=np.radians(np.arange(0,360,1.0));N=np.outer(np.cos(c),e1)+np.outer(np.sin(c),e2)
def along(v,iso=False):
    Vs=fit.V;fit.V=N
    if iso: T1s=fit.T1;fit.T1=lambda p:np.ones_like(p)
    full=sky(aC,v)  # normalised to circle mean; use full-sky mean instead
    fit.V=Vs;mean=None
    if iso: fit.T1=T1s
    return full
fig,ax=plt.subplots(figsize=(16,9));ang_=np.arange(0,360,1.0);pk=np.degrees(np.arccos(PH@e1))
for (y,lab,col) in ((along(None),f"model only: {r[0]['amp_xPlanck']:.2f}× Planck, hot end at CoE",'C1'),(along(PH,True),f"motion only (370 km/s toward Planck hot): {r[1]['amp_xPlanck']:.2f}× Planck",'C0'),
    (along(PH),f"model + motion toward Planck hot: {r[4]['amp_xPlanck']:.2f}× Planck, hot end {r[4]['offset_from_Planck_hot_deg']:.0f}° from Planck hot",'C3'),(along(vec(*b.x)),f"model + motion chosen to cancel (2-angle fit): {r[5]['amp_xPlanck']:.2f}× Planck",'C2')):
    ax.plot(ang_,(y-y.mean())*T0,col,lw=3,label=lab)
ax.axvline(pk,color='r',lw=2);ax.text(pk+2,ax.get_ylim()[1]*.9,'Planck hot',color='r',fontsize=14);ax.axvline((pk+180)%360,color='b',lw=2);ax.text((pk+180)%360+2,ax.get_ylim()[1]*.9,'Planck cold',color='b',fontsize=14)
ax.axvline(0,color='k',ls='--');ax.text(2,ax.get_ylim()[0]*.9,'CoE',fontsize=14);ax.axvline(180,color='k',ls='--');ax.text(182,ax.get_ylim()[0]*.9,'anti-CoE',fontsize=14)
ax.set_xlim(0,360);ax.set_xticks(range(0,361,30));ax.grid(alpha=.3);ax.tick_params(labelsize=14)
ax.set_xlabel(f'angle along the great circle from CoE (180°, −15°) through Planck hot ({pk:.0f}° from CoE)',fontsize=15);ax.set_ylabel('ΔT (μK), circle mean removed',fontsize=16)
ax.legend(loc='lower center',bbox_to_anchor=(.5,1.0),fontsize=13,frameon=False,title='shrinking recombination radius, with mass, CatWISE slope, mean 1+z ≈ 1100',title_fontsize=14)
fig.savefig('dT_combined.png',dpi=90,bbox_inches='tight')

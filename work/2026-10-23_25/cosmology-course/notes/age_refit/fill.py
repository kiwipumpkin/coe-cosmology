"""INTERPOLATION ASSUMPTION variant: fill the CoE-side gap of the visible ring smoothly, then score. usage: fill.py <json>"""
import os,sys,json,numpy as np,warnings;warnings.filterwarnings('ignore')
fn=sys.argv[1]; d=json.load(open(fn)); law,mode,vm=d['tag'].split('_')
os.environ['R0DENSE']='1'
if vm=='gtc': os.environ['VMODE']='gtc'
sys.argv=['x',law,mode]; exec(open('fit_lc.py').read().split("ck=f'fwd")[0])
from scipy.interpolate import CubicHermiteSpline,CubicSpline
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
pts=[(r['ell'],r.get('z_fwd',r.get('z'))) for r in d['rows'] if (r.get('z_fwd') is not None) or r.get('status')=='rec']
th=np.array([p[0] for p in pts]); S=1/(1+np.array([p[1] for p in pts])); l0=d['l0']
miss=sorted(set(r['ell'] for r in d['rows'])-set(th)); out=dict(tag=d['tag'],src=fn,visible=len(th),gap_deg=[min(miss),max(miss)] if miss else None)
if len(th)<8: print('too few'); sys.exit()
g=np.arange(360.); vis=np.zeros(360,bool)
# periodic unwrap so the gap sits in the middle of the sample run
sh=(th+180-np.mean(miss))%360 if miss else th; o_=np.argsort(sh); xs=sh[o_]; ys=np.log(S[o_])
dy=np.gradient(ys,xs)
def nodip(L):
    t=np.radians(g); X=np.c_[np.ones(360),np.cos(t),np.sin(t)]; c=np.linalg.lstsq(X,L,rcond=None)[0]; return L-X[:,1:]@c[1:]
fig,ax=plt.subplots(2,1,figsize=(16,11))
for meth in ('linear','cubic_hermite','cubic_spline'):
    q=(g+180-np.mean(miss))%360 if miss else g
    if meth=='linear': Lg=np.interp(q,xs,ys)
    elif meth=='cubic_hermite': Lg=CubicHermiteSpline(xs,ys,dy)(q)
    else: Lg=CubicSpline(xs,ys)(q)
    inside=(q>=xs[0])&(q<=xs[-1]); Lg=np.where(inside,Lg,np.interp(q,xs,ys))
    res={}
    for t_,LL in (('raw',Lg),('nodip',nodip(Lg))):
        SS=np.exp(LL); So=np.interp((g-l0)%360,np.arange(361.),np.append(SS,SS[0])); rr=CC.resid(CC.G_FULL,So); v=[score(SS,p)[0] for p in np.arange(0,360,1.)]
        res[t_]=dict(cancel_l0=100*(1-CC.rms(rr)/S0),res_l0_1e5=CC.rms(rr)/1e-5,cancel_free=100*(1-min(v)/S0),p2p_dTT=float(np.ptp(So/So.mean())),
            modes={int(m['m']):float(m['A']) for m in CC.fourier(So/So.mean()-1,6)['modes']})
        ax[0 if t_=='raw' else 1].plot(g,So/So.mean()-1,lw=3,label=f"{meth} ({res[t_]['cancel_l0']:.3g}%)")
    out[meth]=res
for a,t_ in zip(ax,('raw','m=1 removed')):
    a.plot(g,CC.G_FULL-1,'k',lw=1.5,label='G−1 (target)'); a.set_ylabel(f'dT/T {t_}',fontsize=18); a.legend(fontsize=14); a.tick_params(labelsize=14); a.grid(); a.set_xlim(0,360)
    if miss: [a.axvspan((m+l0)%360-3,(m+l0)%360+3,color='0.85') for m in miss]
ax[1].set_xlabel('galactic longitude ℓ (deg); grey = filled gap',fontsize=18)
fig.suptitle(f"GAP FILL (interpolation assumption) {d['tag']} [{os.path.basename(fn)}]",fontsize=18); fig.tight_layout()
b=os.path.basename(fn).replace('.json',''); fig.savefig(f'gr2/fill_{b}.png',dpi=80)
json.dump(out,open(f'gr2/fill_{b}.json','w'),indent=1)
print(b,{m:(round(out[m]['raw']['cancel_l0'],1),round(out[m]['nodip']['cancel_l0'],1),'%.3g'%out[m]['nodip']['p2p_dTT']) for m in ('linear','cubic_hermite','cubic_spline')})

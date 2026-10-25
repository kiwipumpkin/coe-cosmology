"""Refit eta fwd gtc Migkas chi2 over all other params (u,p,lf,lL,s) with lM held fixed. usage: refitM.py lM"""
import os,sys,json,numpy as np,warnings;warnings.filterwarnings('ignore')
lm=float(sys.argv[1]); os.environ['R0DENSE']='1'; os.environ['VMODE']='gtc'; sys.argv=['x','eta','fwd']
exec(open('fit_lc.py').read().split("ck=f'fwd")[0]); import grltb as G
from scipy.optimize import minimize
x0=np.array(json.load(open('gr2/v4/best/best_eta_fwd_gtc_1.json'))['x'])
def c(y):
    z=np.r_[lm,y]
    try:
        P,o=build(z)
        if o is None: return 1e3
        f=G.fields(o); af=G.migkas(G.trace(f,np.arange(0,360,4.),nstep=320,lbmax_frac=0.999))
        return F.chi2_l0(af)[0] if np.all(np.isfinite(af)) else 1e3
    except (Exception,MemoryError): return 1e3
# start: find a feasible u by scanning (others held), then polish all five
best=(1e3,x0[1:])
for du in np.linspace(-3,3,25):
    y=x0[1:].copy(); y[0]+=du; v=c(y)
    if v<best[0]: best=(v,y)
r=minimize(c,best[1],method='Nelder-Mead',options=dict(maxiter=150,adaptive=True))
y=r.x if r.fun<best[0] else best[1]; f=min(r.fun,best[0])
json.dump(dict(lM=lm,x=list(np.r_[lm,y]),chi2=float(f)),open(f'gr2/v4/shadow/refit_lM{lm:.2f}.json','w')); print('REFIT',lm,f,list(y))

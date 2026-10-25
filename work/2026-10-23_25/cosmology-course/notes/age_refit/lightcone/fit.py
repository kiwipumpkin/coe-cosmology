"""usage: fit.py <law>.  DE (search box only) -> unbounded Nelder-Mead polish; 12-deg grid for search, 6-deg for every reported number."""
import sys,json,os,numpy as np; from lc import *
from scipy.optimize import differential_evolution,minimize
law=sys.argv[1]; nm=NAMES[law]
BOX={'lM':(-4,1),'lrho':(-4,1),'lL':(-2,2),'s':(0.3,5),'lV':(0,1.5),'p':(-1.5,2),'ltn':(-1.5,1)}
def ev(z,ells=np.arange(0,360,12.)):
    try:
        L,tn=make(law,z); o=obs(L,tn)
        if not np.isfinite(o['U']) or o['U']<=1: return None
        af=migkas(L,tn,ells)
        if af is None or not np.all(np.isfinite(af)): return None
        return (*F.chi2_l0(af),af)
    except Exception: return None
import time
_n=[0]
def cost(z):
    t0=time.time(); r=ev(z); _n[0]+=1; c=1e3 if r is None else r[0]
    print(_n[0],'%.1fs'%(time.time()-t0),round(c,3),flush=True); return c
ck=f'ck_{law}.json'
if __name__=='__main__':
    st=json.load(open(ck)) if os.path.exists(ck) else {}
    if 'de' not in st:
        r=differential_evolution(cost,[BOX[n] for n in nm],seed=3,popsize=6,maxiter=12,polish=False,init='sobol',workers=1,
            callback=lambda xk,convergence=0: json.dump(dict(partial=list(map(float,xk))),open(ck+'.part','w')))
        st['de']=dict(x=list(map(float,r.x)),fun=float(r.fun)); json.dump(st,open(ck,'w')); print('DE',st['de'],flush=True)
    if 'nm' not in st:
        r=minimize(cost,st['de']['x'],method='Nelder-Mead',options=dict(maxiter=60,adaptive=True))
        x=r.x if r.fun<st['de']['fun'] else st['de']['x']; st['nm']=dict(x=list(map(float,x))); json.dump(st,open(ck,'w'))
    x=np.array(st['nm']['x']); chi,l0,af=ev(x,np.arange(0,360,6.))
    st['final6']=dict(chi2=chi,l0=l0,af=af.tolist(),sigma=sig(chi,len(nm)+1)); json.dump(st,open(ck,'w')); print('FINAL',law,chi,l0,flush=True)

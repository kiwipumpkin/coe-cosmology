"""usage: fit2.py <law>: Migkas fit with free lRi (initial displacement), seeded from res_<law>.json; checkpoints."""
import sys,json,os,time,numpy as np; from lc2 import *
from scipy.optimize import minimize
law=sys.argv[1]; nm=NAMES[law]+['lRi']; r0=json.load(open(f'res_{law}.json')); x0=np.array(r0['x']+[-3.0])
def ev(z,ells=np.arange(0,360,12.)):
    try:
        L,tn=make2(law,z[:-1],10**z[-1]); o=obs(L,tn)
        if not np.isfinite(o['U']) or o['U']<=1: return None
        af=migkas(L,tn,ells); return None if af is None or not np.all(np.isfinite(af)) else (*F.chi2_l0(af),af)
    except Exception: return None
n=[0]
def cost(z):
    r=ev(z); c=1e3 if r is None else r[0]; n[0]+=1; print(n[0],round(c,3),list(np.round(z,3)),flush=True)
    json.dump(dict(last=list(map(float,z))),open(f'ck2_{law}.part','w')); return c
ck=f'ck2_{law}.json'; st=json.load(open(ck)) if os.path.exists(ck) else {}
if 'nm' not in st:
    s=np.r_[np.full(len(x0)-1,0.05),1.0]
    rr=minimize(cost,x0,method='Nelder-Mead',options=dict(maxiter=80,adaptive=True,initial_simplex=np.vstack([x0]+[x0+np.eye(len(x0))[i]*s[i] for i in range(len(x0))])))
    st['nm']=list(map(float,rr.x)); json.dump(st,open(ck,'w'))
x=np.array(st['nm']); chi,l0,af=ev(x,np.arange(0,360,6.))
out=dict(law=law,names=nm,x=list(x[:-1]),lRi=float(x[-1]),chi2=chi,l0=l0,sigma=sig(chi,len(nm)+1),af=list(af))
L,tn=make2(law,x[:-1],10**x[-1]); o=obs(L,tn); lam=o['Hbar']*F.C/(70/F.KMSMPC)
out.update(U_c=o['U'],D_Gpc=o['R']*lam,Hpar=70*o['Hpar']/o['Hbar'],Hper=70*o['Hper']/o['Hbar'],age_Gyr=tn/F.C*lam,M_kg=L.P['M']*lam)
json.dump(out,open(f'res2_{law}.json','w'),indent=1,default=float); print('FINAL',chi,out['sigma'],flush=True)

# Phase5 (NEW): Migkas A/A_all (A=af^2/<af^2>) refit at FIXED d. Free: D, mc, l0, eta, k (L derived from CatWISE lock delta=0.031k).
# <H0>=68.64 locked; age 13.787 Gyr anchor (5% penalty, same as phase4b). No parameter box: unbounded Nelder-Mead in log params;
# only physical guards (k4b.sky: d<D/1.05, L root exists, vrmin>0). Multistart (phase4b best, phase4 best, 4 random perturbations).
# Usage: python3 fitd5.py d  -> fitd_<d>.json (per-start checkpoint fitd_<d>_ckpt.json)
import sys,os,json,numpy as np,warnings;warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/coepdf/joint_solve/phase4b');import k4b as K
from scipy.optimize import minimize
O='/workspace/coepdf/joint_solve/phase5/';d=float(sys.argv[1]);TOL=0.05
def parts(x):
    D,mc,l0,eta,k=10**x[0],10**x[1],x[2],10**x[3],10**x[4];o=K.sky(D,d,eta,mc,k)
    if o is None:return None
    a2=o[0]**2/np.mean(o[0]**2);m=K.cw.shifted(a2,l0%360)[K.cw.IMIG];return float(np.sum(K.cw.W*(m-K.cw.Y)**2)),o[1]['Ht']*K.TH,m,o
def obj(x):
    try:p=parts(x)
    except Exception:return 800.
    return 800. if p is None else p[0]+((p[1]/K.AGE-1)/TOL)**2
S=[[1.4368,-1.252,18.64,1.0607,3.1002],[1.4459,-1.1866,18.81,0.3422,1.9997]]
rng=np.random.default_rng(int(d*1000)+1)
for i in range(4):S.append((np.array(S[i%2])+rng.normal(0,[0.03,0.3,2,0.4,1])).tolist())
ck=O+f'fitd_{d:g}_ckpt.json';C=json.load(open(ck)) if os.path.exists(ck) else []
for s in S[len(C):]:
    r=minimize(obj,np.array(s),method='Nelder-Mead',options=dict(maxfev=int(os.environ.get('FEV','1500')),xatol=1e-5,fatol=1e-5,adaptive=True))
    C.append(dict(x0=s,x=r.x.tolist(),f=float(r.fun)));json.dump(C,open(ck,'w'));print(d,'start',round(r.fun,4),flush=True)
b=min(C,key=lambda c:c['f']);x=np.array(b['x'])
r=minimize(obj,x,method='Nelder-Mead',options=dict(maxfev=1500,xatol=1e-6,fatol=1e-6,adaptive=True));x=r.x if r.fun<b['f'] else x
chi,age,m,o=parts(x);Y=K.cw.Y
R=dict(d=d,chi2=chi,obj=obj(x),r2=float(1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)),D=10**x[0],mc=10**x[1],l0=x[2]%360,eta=10**x[3],k=10**x[4],delta=0.031*10**x[4],L=o[2]*10**x[0],age_Gyr=age,Ht=o[1]['Ht'],n=int(len(Y)),
 starts=[c['f'] for c in C],flags=[s for s,v in [('eta>500',10**x[3]>500),('k>5e4',10**x[4]>5e4),('k<1e-2',10**x[4]<1e-2),('mc<1e-5',10**x[1]<1e-5),('D near 1.05d guard',10**x[0]<1.1*d)] if v])
json.dump(R,open(O+f'fitd_{d:g}.json','w'),indent=1);print('DONE',R,flush=True)

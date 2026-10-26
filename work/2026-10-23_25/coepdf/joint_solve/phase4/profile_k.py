# Phase4 (NEW): profile likelihood of Migkas chi2 (+age penalty, TOL=0.05) vs gradient ratio k = delta/0.031.
# Two versions: (A) CatWISE depth d free in [2,40] Gpc; (B) d fixed at the code's CatWISE lock depth RV=14.30 Gpc.
import numpy as np,json,os,sys,multiprocessing as mp
from scipy.optimize import minimize
sys.path.insert(0,'/workspace/coepdf/joint_solve/phase4')
def job(a):
    import joint_fit as J
    mode,lk,x0=a;x0=np.array(x0);dfix=np.log10(J.K.cw.RV)
    J.B=J.B0[:5]+[(-2,3)]
    def f(y):
        x=np.r_[y[0],dfix,y[1:],lk] if mode=='B' else np.r_[y,lk];return J.obj(x)
    y0=np.r_[x0[0],x0[2:5]] if mode=='B' else x0[:5];best=None
    for s in (y0,y0+np.r_[.1,0,0,.3][:len(y0)] if mode=='B' else y0+np.r_[.1,.2,0,0,.3]):
        r=minimize(f,s,method='Nelder-Mead',options=dict(maxfev=800,xatol=1e-4,fatol=1e-4,adaptive=True))
        if best is None or r.fun<best.fun:best=r
    x=np.r_[best.x[0],dfix,best.x[1:],lk] if mode=='B' else np.r_[best.x,lk]
    try:R=J.rep(x)
    except Exception:R=dict(fail=True)
    return mode,lk,float(best.fun),{k:float(v) for k,v in R.items()}
if __name__=='__main__':
    bx=json.load(open('/workspace/coepdf/joint_solve/phase4/joint_T0.05.json'))['best_x']
    jobs=[(m,lk,bx) for m in 'AB' for lk in np.linspace(-2,3,21)]
    out=[]
    with mp.get_context('spawn').Pool(8) as p:
        for r in p.imap_unordered(job,jobs):out.append(r);print(r[0],round(r[1],2),round(r[2],3),flush=True)
    json.dump(sorted(out),open('/workspace/coepdf/joint_solve/phase4/profile_k.json','w'),indent=1)

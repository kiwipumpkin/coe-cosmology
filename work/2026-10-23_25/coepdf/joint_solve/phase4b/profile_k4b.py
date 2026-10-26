# Phase4b (NEW): profile chi2(+age penalty) vs k=delta/0.031, all other params free in joint_fit4b bounds; (A) d free, (B) d=14.30 Gpc.
import numpy as np,json,sys,multiprocessing as mp
from scipy.optimize import minimize
O='/workspace/coepdf/joint_solve/phase4b/'
def job(a):
    sys.path.insert(0,O);import joint_fit4b as J
    mode,lk,x0=a;x0=np.array(x0);df=np.log10(14.30048287260589)
    def full(y):return np.r_[y[0],df,y[1:],lk] if mode=='B' else np.r_[y,lk]
    y0=np.r_[x0[0],x0[2:5]] if mode=='B' else x0[:5];best=None
    seeds=[y0]+[y0+np.r_[e] for e in ([.1,.3,0,.3],[-.1,-.3,0,-.3],[.3,.6,0,.6]) ] if mode=='B' else [y0]+[y0+np.r_[e] for e in ([.1,.3,.3,0,.3],[-.1,-.5,-.3,0,-.3],[.3,-1,.6,0,.6],[0,.5,0,0,0])]
    for s in seeds:
        r=minimize(lambda y:J.obj(full(y)),s,method='Nelder-Mead',options=dict(maxfev=900,xatol=1e-4,fatol=1e-4,adaptive=True))
        if best is None or r.fun<best.fun:best=r
    try:R=J.rep(full(best.x))
    except Exception:R=dict(fail=1.)
    return mode,float(lk),float(best.fun),{k:float(v) for k,v in R.items()}
if __name__=='__main__':
    bx=json.load(open(O+'joint4b.json'))['best_x'];jobs=[(m,lk,bx) for m in 'AB' for lk in np.linspace(-1,5,25)];out=[]
    with mp.get_context('spawn').Pool(8) as p:
        for r in p.imap_unordered(job,jobs):out.append(r);print(r[0],round(r[1],2),round(r[2],3),flush=True);json.dump(sorted(out),open(O+'profile_k4b.json','w'),indent=1)

import sys,pickle,warnings; warnings.filterwarnings("ignore"); sys.path.insert(0,'.'); sys.path.insert(0,'..')
import run_exactc as X, ltbcore, numpy as np, multiprocessing as mp
from ltbcore import LTB, worldline, photon
from grcore import stats, norm
import functools
def ok180(L,v0,T):
    g=1/np.sqrt(1-v0*v0); oat,_=worldline(L,X.C0+np.array([g*v0*v0,0]),np.array([v0,0]),L.t1+g*v0,T); y=oat(T)
    so=photon(L,y[0],y[1],0.0,T,L.t1-0.1)
    if so.t[-1]>L.t1: return False
    return so.sol(L.t1)[0] < X.D-1/g
def job(a):
    Rdot,core,v0=a
    try:
        orig=X.LTB; X.LTB=functools.partial(orig,Rdot=Rdot)
        L=X.LTB(core_frac=0.0946,rmax=150.0,nr=900); lo=L.t1+2*v0/np.sqrt(1-v0*v0)+0.05; hi=lo+2
        for _ in range(22):
            m=(lo+hi)/2; (hi:=m) if ok180(L,v0,m) else (lo:=m)
        dt=np.ceil((hi-L.t1)/0.1-1e-9)*0.1
        for k in range(4):
            x=X.run(core,v0,L.t1+dt)
            if np.isnan(x).sum()<=6: return a,(dt,x,L.t1)
            dt+=0.1
        return a,(dt,x,L.t1)
    except Exception as ex: return a,repr(ex)
if __name__=='__main__':
    jobs=[(R,c,v) for R in (1.5,6.0) for v in (0.0,0.00123,0.1) for c in (0.0946,0.0)]
    out={}
    with mp.Pool(6) as P:
        for a,r in P.imap_unordered(job,jobs):
            out[a]=r; pickle.dump(out,open('recession.pkl','wb'))
            print(a, r if isinstance(r,str) else (r[0], r[2], int(np.isnan(r[1]).sum()), np.round(np.array(stats(norm(np.nan_to_num(r[1],nan=np.nanmean(r[1])))))[[0,1,3]],5)),flush=True)

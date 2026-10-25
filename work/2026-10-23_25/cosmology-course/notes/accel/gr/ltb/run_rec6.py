import sys,pickle; import run_recession as RR, run_exactc as X, functools, numpy as np, multiprocessing as mp
from grcore import stats, norm
F=0.0946/4
def job(a):
    core,v0=a; X.LTB=functools.partial(RR.LTB,Rdot=6.0,rmin=7.9)
    try:
        L=X.LTB(core_frac=F,rmax=150.0,nr=900); lo=L.t1+2*v0/np.sqrt(1-v0*v0)+0.05; hi=lo+2
        for _ in range(22):
            m=(lo+hi)/2; (hi:=m) if RR.ok180(L,v0,m) else (lo:=m)
        dt=np.ceil((hi-L.t1)/0.1-1e-9)*0.1; x=X.run(core,v0,L.t1+dt); return a,(dt,x,L.t1)
    except Exception as ex: return a,repr(ex)
if __name__=='__main__':
    out={}
    with mp.Pool(3) as P:
        for a,r in P.imap_unordered(job,[(F,0.0),(F,0.00123),(F,0.1)]):
            out[a]=r; pickle.dump(out,open('recession6_fixedmass.pkl','wb'))
            print(a, r if isinstance(r,str) else (r[0],int(np.isnan(r[1]).sum()),np.round(np.array(stats(norm(np.nan_to_num(r[1],nan=np.nanmean(r[1])))))[[0,1,3]],5)),flush=True)

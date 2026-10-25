# Point-mass LTB CMB ring: exactly-c shell => H*RV=c, so D/RV = space speed at D (in c). CoE at ring angle 180.
import sys,os,pickle,functools,warnings;warnings.filterwarnings('ignore')
G='/workspace/cosmology-course/notes/accel/gr';sys.path[:0]=[G+'/valley',G+'/ltb',G]
os.chdir(G+'/valley');import toy;os.chdir(os.path.dirname(os.path.abspath(__file__)))
import numpy as np, run_exactc as X
F=0.0178
def job(a):
    uD,f=a; Dt=uD; Mc=max(f,1e-12)*Dt**3/2
    Mfun=lambda r: Mc+0*r
    rm=max(0.02,Dt-1.0-0.03)
    while True:
        try: L=toy.LTBM(Mfun,Dt,1.0*Dt,rm,rmax=150.0,nr=900); break
        except ValueError: rm+=0.02
        if rm>Dt-0.5: return a,'no grid'
    X.LTB=lambda core_frac=None,rmax=150.0,nr=900: toy.LTBM(Mfun,Dt,1.0*Dt,rm,rmax=rmax,nr=nr)
    X.D=Dt;X.C0=np.array([Dt,0.0])
    for dt in (1.3,1.8,2.5,3.5):
        try: x=X.run(None,0.00123,L.t1+dt)
        except Exception as e: x=None;err=repr(e);continue
        if not np.isnan(x).any(): break
    if x is None: return a,err
    return a,dict(x=x,dt=dt,rmin=rm,t1=L.t1,twoM_over_D=2*Mc/Dt)
if __name__=='__main__':
    import multiprocessing as mp
    jobs=[(u,f) for u in (1.05,3.0,10.5) for f in (F,0.0)]
    out={}
    with mp.Pool(6) as P:
        for a,r in P.imap_unordered(job,jobs):
            out[a]=r;pickle.dump(out,open('pcmb.pkl','wb'))
            print(a,r if isinstance(r,str) else (r['dt'],r['rmin'],int(np.isnan(r['x']).sum()),np.nanstd(r['x'])),flush=True)

# Valley point -> LTB CMB toy (exactly-c shell, observer at 0.00123c). Ring radius = vision sphere RV; length unit RV; H*RV/c = 0.3 (assumption)
import sys,json,pickle,functools,warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/ltb'); sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr')
import numpy as np, ltbcore
from ltbcore import *
import run_exactc as X, run_recession as RR
BUB="/workspace/cosmology/papers/spreadsheet_model/bubble"; sys.path.insert(0,BUB); import os; _c=os.getcwd(); os.chdir(BUB); import ptshell, prof_core as P; os.chdir(_c)
RV=14.30048287260589; OS=0.049; HR=0.3
class LTBM(LTB):
    def __init__(s,Mfun,D,Rdot,rmin,rmax=150.0,tspan=400.0,nr=900,nt=700,**kw):
        s.D=D; M=Mfun; s.Mf=M; s.Mc=0.0
        ED=(Rdot**2-2*M(D)/D)/2; s.t1=s.age(D,M(D),ED)
        rg=np.linspace(rmin,rmax,nr); Eg=np.array([s.solveE(r,M(r)) for r in rg]); s.Es=CubicSpline(rg,Eg)
        tg=s.t1+np.r_[-s.t1*(1-np.geomspace(0.0005,1,160))[:-1], np.geomspace(1e-3,tspan,nt-6)]
        tg=np.unique(tg); tg=tg[tg>0.0004*s.t1]; Rg=np.empty((nr,len(tg)))
        for i,(r,E) in enumerate(zip(rg,Eg)):
            Mr=M(r); f=lambda t,R:[np.sqrt(max(2*Mr/R[0]+2*E,1e-30))]
            so=solve_ivp(f,(s.t1,tg[-1]),[r],t_eval=tg[tg>=s.t1],rtol=1e-12,atol=1e-12)
            sb=solve_ivp(f,(s.t1,tg[0]),[r],t_eval=tg[tg<s.t1][::-1],rtol=1e-12,atol=1e-12)
            Rg[i]=np.r_[sb.y[0][::-1],so.y[0]]
        s.Rs=RectBivariateSpline(rg,tg,Rg,kx=5,ky=5); s.rg=rg; s.tg=tg; s.Eg=Eg
        s.Rr=s.Rs.partial_derivative(1,0); s.Rt_=s.Rs.partial_derivative(0,1); s.Rrt=s.Rs.partial_derivative(1,1); s.Rrr=s.Rs.partial_derivative(2,0)
        s.minRp=s.Rs(rg,tg,dx=1).min(); s.min1p2E=(1+2*Eg).min()
def make(p,core=True):
    Dt=p['D']/RV; a=dict(DR=Dt,Os=OS,mcc=p['mcc'] if core else 1e-12,rin=1-p['kin']/Dt,rout=1+p['kout']/Dt,Lc=0.005)
    S,MT,m1=ptshell.build(a); Mcrit=HR**2*Dt**3/2
    Mfun=lambda r: m1*Mcrit*np.interp(r/Dt,P.RG,MT)
    rmin=max(0.15,Dt-2.6)
    return functools.partial(LTBM,Mfun,Dt,HR*Dt,rmin), Dt
def job(i):
    p=json.load(open('valley_points.json'))[i]; v0=0.00123
    try:
        fac,Dt=make(p); rm=max(0.15,Dt-2.6)
        while True:
            try: L=fac.func(*fac.args[:3],rm); break
            except ValueError: rm+=0.05
        X.LTB=lambda core_frac=None,rmax=150.0,nr=900: fac.func(*fac.args[:3],rm,rmax=rmax,nr=nr)
        X.D=Dt; X.C0=np.array([Dt,0.0])
        L=X.LTB()
        for dt in (1.3,1.8,2.5):
            x=X.run(None,v0,L.t1+dt)
            if not np.isnan(x).any(): break
        return i,dict(rmin=rm,dt=dt,t1=L.t1,minRp=L.minRp,min1p2E=L.min1p2E,x=x)
    except Exception as ex: return i,repr(ex)
if __name__=='__main__':
    import multiprocessing as mp
    from grcore import stats, norm
    idx=[int(a) for a in sys.argv[1:]]; out={}
    with mp.Pool(len(idx)) as Pp:
        for i,r in Pp.imap_unordered(job,idx):
            out[i]=r; pickle.dump(out,open('toy_%s.pkl'%'_'.join(sys.argv[1:]),'wb'))
            print(i, r if isinstance(r,str) else (r['dt'],round(r['minRp'],3),round(r['min1p2E'],3),int(np.isnan(r['x']).sum()),np.round(np.array(stats(norm(np.nan_to_num(r['x'],nan=np.nanmean(r['x'])))))[[0,1,3]],5)),flush=True)

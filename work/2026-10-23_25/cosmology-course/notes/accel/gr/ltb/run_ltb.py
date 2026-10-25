import sys,pickle,warnings; warnings.filterwarnings("ignore"); sys.path.insert(0,'.'); sys.path.insert(0,'..')
from ltbcore import *; import multiprocessing as mp
from grcore import relsum, stats, norm
th=np.radians(np.arange(0,361,1.0)); D=10.0; C0=np.array([D,0.0]); U=0.9999
LL={}
def getL(core):
    if core not in LL: LL[core]=LTB(core_frac=core,rmax=150.0,nr=900)
    return LL[core]
def setup(L,v0,T):
    g=1/np.sqrt(1-v0*v0); ems=[]; st=[]
    for t in th:
        e=np.array([np.cos(t),np.sin(t)]); ts=L.t1+g*v0*(1+e[0])
        at,s_=worldline(L,C0+np.array([g*(e[0]+v0*v0),e[1]]),relsum(v0,U,e),ts,T); ems.append((at,ts)); st.append(s_)
    oat,os_=worldline(L,C0+np.array([g*v0*v0,0]),np.array([v0,0]),L.t1+g*v0,T)
    return ems,oat,sum(1 for x in st if x==1),os_
def ok180(L,v0,T):
    g=1/np.sqrt(1-v0*v0); oat,_=worldline(L,C0+np.array([g*v0*v0,0]),np.array([v0,0]),L.t1+g*v0,T); y=oat(T)
    if y is None: raise RuntimeError('observer left the grid')
    so=photon(L,y[0],y[1],0.0,T,L.t1-0.1)
    if so.t[-1]>L.t1: return False
    return so.sol(L.t1)[0] < D-1/g
def run(core,v0,T):
    L=getL(core); ems,oat,nedge,oedge=setup(L,v0,T); y=oat(T)
    if y is None: return None,nedge,"observer left the grid"
    return observe(L,y,ems,T),nedge,None
def job(v0):
  try: return job_(v0)
  except Exception as ex: return v0,None,None,repr(ex)
def job_(v0):
    L=getL(0.0946); g=1/np.sqrt(1-v0*v0); lo=L.t1+2*g*v0+0.05; hi=lo+2
    while not ok180(L,v0,hi):
        hi=L.t1+(hi-L.t1)*1.5
        if hi-L.t1>350: return v0,None,None,"light from behind had not caught up by 350 R/c"
    for _ in range(22):
        mid=(lo+hi)/2
        if ok180(L,v0,mid): hi=mid
        else: lo=mid
    dt=hi-L.t1; step=0.1 if dt<5 else (0.5 if dt<30 else 1.0); dt=np.ceil(dt/step-1e-9)*step
    for k in range(15):
        res={c:run(c,v0,L.t1+dt) for c in (0.0,0.0946)}
        if all(r[0] is not None and not np.isnan(r[0]).any() for r in res.values()): return v0,dt,res,None
        dt+=step
    return v0,dt,res,"still missing directions"
if __name__=='__main__':
    V=[float(a) for a in sys.argv[1:]]
    with mp.Pool(min(6,len(V))) as P:
        out={}
        for v0,dt,res,flag in P.imap_unordered(job,V):
            out[v0]=(dt,res,flag); pickle.dump(out,open('ltb_%s.pkl'%'_'.join(sys.argv[1:]),'wb'))
            def s(r):
                if r[0] is None: return r[2]
                x=r[0]
                if np.isnan(x).any(): return "missing %d"%np.isnan(x).sum()
                d,q,_,mn=stats(norm(x)); return "%.4f %.4f %.4f (edge %d)"%(d,q,mn,r[1])
            print(v0,"dt",dt,flag,"" if res is None else " | ".join(f"core {c}: "+s(r) for c,r in res.items()),flush=True)

# Shell edge receding at EXACTLY c: emitters are null geodesics (U->1 limit of the old setup, gamma factored out)
import sys,pickle,warnings; warnings.filterwarnings("ignore"); sys.path.insert(0,'.'); sys.path.insert(0,'..')
from ltbcore import *; import ltbcore; import multiprocessing as mp
from grcore import relsum, stats, norm
th=np.radians(np.arange(0,361,1.0)); D=10.0; C0=np.array([D,0.0])
def kdot(L,t,r,vk,lnk,vl,lnl):   # k.l for two null vectors given coordinate velocities and ln(time components)
    R,Rp,Rt,G,_,_=L.fields(t,r); return np.exp(lnk+lnl)*(-1+G*vk[0]*vl[0]+R*R*vk[1]*vl[1])
def run(core,v0,T):
    L=LTB(core_frac=core,rmax=150.0,nr=900); g=1/np.sqrt(1-v0*v0); ems=[]; l0=[]
    for t in th:
        e=np.array([np.cos(t),np.sin(t)]); ts=L.t1+g*v0*(1+e[0])
        w=relsum(v0,1.0,e); w/=np.linalg.norm(w)          # exactly light-like, in the local comoving frame
        at,_=worldline(L,C0+np.array([g*(e[0]+v0*v0),e[1]]),w,ts,T); ems.append((at,ts))
        l0.append(g*(1+v0*e[0]))   # time component of the emitter's null vector in the comoving frame at start = limit of u^t/gamma_U
    oat,_=worldline(L,C0+np.array([g*v0*v0,0]),np.array([v0,0]),L.t1+g*v0,T); y=oat(T)
    # same ray matching as ltbcore.observe, but emitter factor (k.l)_e instead of (k.u)_e
    ltbcore_kdu=L.kdu; out=[]
    for i,(em,ts) in enumerate(ems):
        idx=[i]
        def fake_kdu(t,r,vph,lnk,vp,_i=i,_em=em):
            if t==T and lnk==0.0: return ltbcore_kdu(t,r,vph,lnk,vp)
            e=_em(t); return kdot(L,t,r,vph,lnk,e[2:4],e[4]+np.log(l0[_i]))
        L.kdu=fake_kdu; out.append(observe(L,y,[(em,ts)],T)[0]); L.kdu=ltbcore_kdu
    return np.array(out)
def job(a):
    core,v0,dt=a
    try: return a,run(core,v0,2.2222222222222223+dt)
    except Exception as ex: return a,repr(ex)
if __name__=='__main__':
    DT={0.0:1.4,0.00123:1.4,0.01:1.4,0.1:1.6,0.5:2.4}
    jobs=[(c,v,DT[v]) for v in DT for c in (0.0946,0.0)]
    out={}
    with mp.Pool(6) as P:
        for a,x in P.imap_unordered(job,jobs):
            out[a]=x; pickle.dump(out,open('exactc.pkl','wb'))
            s=x if isinstance(x,str) else ("missing %d"%np.isnan(x).sum() if np.isnan(x).any() else "%.5f %.5f %.5f"%tuple(np.array(stats(norm(x)))[[0,1,3]]))
            print(a,s,flush=True)

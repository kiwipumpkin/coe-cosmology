"""Bodies on top of a FIXED background, original joint approach (asym_core-style, no-mirror): all N bodies fitted together
(mass or void radius, l, b, distance) by bounded least squares from multiple random starts. No greedy, no candidate pools.
Target: T = R (SR + central-mass residual) - w (the background's ring pattern).
Region (standard): our sphere + shell out to 30 Gpc from us, including the CoE-core cap inside 30 Gpc. Bounded positions:
  distance from us d in [R_body + 0.5, 30 - R_body]; any direction.
Hard rules (as penalties driven to zero, then checked): every body's SURFACE >= 3 Gpc from the CMB ring (centre distance from
  the ring circle minus body radius >= 3); no straddling; bodies don't overlap; no-mirror rule of asym_core (no same-kind pair
  within 20 deg / 15 % distance / 25 % size of each other's mirror image about the CoE plane; nothing within 10 deg of that plane).
Bodies: clusters 1e49-2e52 kg (compact, delta_c), deep voids delta -0.7, R 0.5-4 Gpc. Compensation: all-open FLOC convention.
Starts: positions uniform in VOLUME over the allowed region, sizes uniform in the size parameter.
usage: python3 joint_add.py NPROC  (job list built in-file; results -> joint_add/raw/<tag>_N<N>_s<start>.json, skipped if present)"""
import os,sys,json,time,numpy as np
os.environ.setdefault("OMP_NUM_THREADS","1")
from multiprocessing import Pool
from scipy.optimize import least_squares
WEBDIR=os.path.dirname(os.path.abspath(__file__)); CO="/workspace/clump_scan/core_open"; sys.path.insert(0,CO)
_c=os.getcwd(); os.chdir(CO); from co_core import *; import ao_core; os.chdir(_c)
from fit_on_web import spec
RT=30.0; CLR=float(os.environ.get("RING_MIN",3.0)); ANG=20.0; DT=np.log(1.15); STl=np.log(1.25); PLANE=10.0
AXR=np.radians(15.16); A2=2*AXR; CA,SA=np.cos(A2),np.sin(A2); NRM=np.array([-np.sin(AXR),np.cos(AXR),0.0])
SC=1/(FLOOR*np.sqrt(360))
def dring(X): return np.sqrt((np.hypot(X[:,0],X[:,1])-RMIN)**2+X[:,2]**2)
class J(ao_core.AO):
    def decode(self,p):
        p=p.reshape(self.n,4); u,lon,lat,s=p.T; isC=self.k=="C"
        M=np.zeros(self.n); Rv=np.zeros(self.n); Rc=np.zeros(self.n)
        M[isC]=10**(np.log10(MLO)+s[isC]*np.log10(MHI/MLO)); Rv[~isC]=RLO+s[~isC]*(RHI-RLO); Rc[~isC]=Rv[~isC]
        Rc[isC]=(3*M[isC]/(4*np.pi*DELTA_C*RHO*GPC**3))**(1/3)
        dmin=Rc+0.5; dmx=RT-Rc; d=(dmin**3+u*(dmx**3-dmin**3))**(1/3)        # u uniform -> uniform in volume
        X=np.c_[np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)]*d[:,None]
        isV=self.k=="V"; M[isV]=-0.70*rho_gpc(X[isV])*vol(Rv[isV])*GPC**3
        return X,M,Rv,Rc
    def violations(self,p):
        X,M,Rv,Rc=self.decode(p); v=[super().violations(p), np.maximum(0,CLR-(dring(X)-Rc))]
        v.append(self.mirror(X,M,Rv)); return np.concatenate(v)
    def mirror(self,X,M,Rv):
        d=np.linalg.norm(X,axis=1); U=X/d[:,None]; Ua=U.copy(); Ua[:,2]=np.abs(Ua[:,2])
        Mi=np.c_[CA*Ua[:,0]+SA*Ua[:,1],SA*Ua[:,0]-CA*Ua[:,1],Ua[:,2]]
        size=np.where(self.k=="C",np.abs(M),Rv); i,j=self.pi,self.pj; out=[]
        if len(i):
            ch=np.linalg.norm(Mi[i]-Ua[j],axis=1); ang=np.degrees(2*np.arcsin(np.clip(ch/2,0,1)))
            s=((ang/ANG)**8+(np.abs(np.log(d[j]/d[i]))/DT)**8+(np.abs(np.log(size[j]/size[i]))/STl)**8)**(1/8); out.append(np.maximum(0,1.17-s))
        plane=np.degrees(np.arcsin(np.clip(np.abs(U@NRM),0,1))); out.append(np.maximum(0,(PLANE*1.05-plane)/PLANE)); return np.concatenate(out)
    def setup(self):
        i,j=np.triu_indices(self.n,1); same=self.k[i]==self.k[j]; self.pi,self.pj=i[same],j[same]; self.lo=np.zeros(4*self.n); self.hi=np.ones(4*self.n)
        self.lo[1::4]=-np.pi; self.hi[1::4]=np.pi; self.lo[2::4]=-np.pi/2; self.hi[2::4]=np.pi/2; return self
def start(m,rng):
    n=m.n; p=np.zeros((n,4)); p[:,3]=rng.uniform(0,1,n)
    for t in range(200):   # uniform-in-volume positions, resampled until the ring rule holds at the start
        p[:,0]=rng.uniform(0,1,n); p[:,1]=rng.uniform(-np.pi,np.pi,n); p[:,2]=np.arcsin(rng.uniform(-1,1,n))
        X,M,Rv,Rc=m.decode(p.ravel()); bad=(dring(X)-Rc)<CLR
        if not bad.any(): break
        p[bad,3]*=0.7
    return p.ravel()
def fit(m,p,T):
    lam=30.0
    for it in range(5):
        f=lambda q: np.concatenate([(T-m.pattern(q))*SC, lam*m.violations(q)])
        r=least_squares(f,p,bounds=(m.lo,m.hi),method="trf",x_scale="jac",max_nfev=40*len(p),xtol=1e-8,ftol=1e-8); p=r.x
        if m.violations(p).max(initial=0)<=1e-3: break
        lam*=10
    return p
def job(a):
    tag,N,s=a; fn=f"{WEBDIR}/joint_add/raw/{tag}_N{N}_s{s}.json"
    if os.path.exists(fn): return fn
    T=R-np.load(f"{WEBDIR}/final_{tag}.npz")["w"]; kinds="C"*(N//2)+"V"*(N-N//2); m=J(kinds).setup(); rng=np.random.default_rng(100*N+s)
    t0=time.time(); p=fit(m,start(m,rng),T); X,M,Rv,Rc=m.decode(p); w=m.pattern(p); res=T-w; v=m.violations(p).max(initial=0)
    out=dict(tag=tag,N=N,start=s,left=float(rms(res)/FLOOR),ok=bool(v<=1e-3),viol=float(v),X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),Rc=Rc.tolist(),
             dist_us=np.linalg.norm(X,axis=1).tolist(),dist_ring=dring(X).tolist(),surface_ring=(dring(X)-Rc).tolist(),in_core=int((core_depth(X)>0).sum()),
             inside_us=int((np.linalg.norm(X,axis=1)<RMIN).sum()),spec=spec(res),background=spec(T),resid=res.tolist(),secs=time.time()-t0)
    json.dump(out,open(fn+".tmp","w")); os.replace(fn+".tmp",fn)
    print(f"{tag} N={N} s={s} left={out['left']:.2f} ok={out['ok']} t={out['secs']:.0f}s",flush=True); return fn
TAGS=["none_s0","shallow_cmb_aware_s2_tidalcore","b_random_s2_tidalcore","b_foam_rim_s3_tidalcore","wave_s1","lattice_s1"]
NS=[8,12,16,20,24,32]; NSTART=int(os.environ.get("NSTART",8))
if __name__=="__main__":
    jobs=[(t,N,s) for t in TAGS for N in NS for s in range(NSTART)]
    # first pass order: none + shallow first, small N first
    jobs.sort(key=lambda j:(TAGS.index(j[0])>1, j[1], TAGS.index(j[0]), j[2]))
    with Pool(int(sys.argv[1])) as P:
        for _ in P.imap_unordered(job,jobs): pass

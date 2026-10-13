"""Fixed background = a filled growth web (random / foam; core / inside / both). Add bodies on top, no restriction on type or place
(masses 1e49..2e52 kg, big or small; deep voids 0.5..4 Gpc; anywhere 0.5..102 Gpc: inside our sphere, shell, core; core_open COAO model),
4 at a time (2 masses + 2 voids), every round re-fitting ALL added bodies jointly, until leftover < 1x or NMAX.
Shallow voids are not offered: in every earlier fit they were inert (no pull outside 3 Rv).
Added bodies don't overlap each other or straddle the ring; they are not checked against the background's spheres.
usage: python3 add_until.py TAG -> web_add/<TAG>.json (resumes from its last round)"""
import os,sys,json,time,numpy as np
os.environ.setdefault("OMP_NUM_THREADS","1")
from scipy.optimize import least_squares
WEBDIR=os.path.dirname(os.path.abspath(__file__)); CO="/workspace/clump_scan/core_open"; sys.path.insert(0,CO)
_c=os.getcwd(); os.chdir(CO); from co_core import *; os.chdir(_c)
sys.path.insert(0,WEBDIR); from fit_on_web import spec
SC=1/(FLOOR*np.sqrt(360)); NMAX=int(os.environ.get("NMAX",120)); TRIALS=3
def fit_new(kold,pold,knew,q0,T):
    """fit only the new block against T - (old bodies), overlap checked against the full set"""
    full=COAO(kold+knew); mnew=COAO(knew); wold=COAO(kold).pattern(pold) if kold else np.zeros(360); lam=30.0; q=q0
    lo=full.lo[len(pold):]; hi=full.hi[len(pold):]
    for it in range(4):
        f=lambda x: np.concatenate([(T-wold-mnew.pattern(x))*SC, lam*full.violations(np.r_[pold,x])])
        r=least_squares(f,q,bounds=(lo,hi),method="trf",x_scale="jac",max_nfev=80*len(q),xtol=1e-9,ftol=1e-9); q=r.x
        if full.violations(np.r_[pold,q]).max(initial=0)<=1e-3: break
        lam*=10
    return np.r_[pold,q]
def refit(mod,p,T):
    lam=30.0
    for it in range(4):
        f=lambda q: np.concatenate([(T-mod.pattern(q))*SC, lam*mod.violations(q)])
        r=least_squares(f,p,bounds=(mod.lo,mod.hi),method="trf",x_scale="jac",max_nfev=40*len(p),xtol=1e-9,ftol=1e-9); p=r.x
        if mod.violations(p).max(initial=0)<=1e-3: break
        lam*=10
    return p
tag=sys.argv[1]; T=R-np.load(f"{WEBDIR}/final_{tag}.npz")["w"]; os.makedirs(f"{WEBDIR}/web_add",exist_ok=True); fn=f"{WEBDIR}/web_add/{tag}.json"
st=json.load(open(fn)) if os.path.exists(fn) else dict(tag=tag,background=spec(T),rounds=[],done=False)
kinds=st["rounds"][-1]["kinds"] if st["rounds"] else ""; p=np.array(st["rounds"][-1]["p"]) if st["rounds"] else np.zeros(0)
rng=np.random.default_rng(abs(hash(tag))%2**32+len(kinds)); t0=time.time()
while not st["done"]:
    newk="CCCCVVVV" if len(kinds)>=8 else "CCVV"; nn=len(newk); k2=kinds+newk; best=None
    for t in range(TRIALS):
        q=np.zeros((nn,4)); q[:,0]=rng.uniform(0,1,nn)**2; q[:,1]=rng.uniform(-np.pi,np.pi,nn); q[:,2]=np.arcsin(rng.uniform(-1,1,nn)); q[:,3]=rng.uniform(0,1,nn)
        mod=COAO(k2); pf=fit_new(kinds,p,newk,q.ravel(),T)
        if len(k2)<=12 or len(k2)%24==0: pf=refit(mod,pf,T)        # joint polish of everything at N<=12 and every 24 bodies
        L=rms(T-mod.pattern(pf))/FLOOR
        if mod.violations(pf).max(initial=0)<=1e-3 and (best is None or L<best[0]): best=(L,pf)
    if best is None: print("no feasible round",flush=True); break
    L,p=best; kinds=k2; mod=COAO(kinds); X,M,Rv,Rc=mod.decode(p); wb=mod.pattern(p); res=T-wb
    st["rounds"].append(dict(N=len(kinds),kinds=kinds,left=float(L),p=p.tolist(),X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),spec=spec(res),resid=res.tolist(),
        inside_us=int((np.linalg.norm(X,axis=1)<RMIN).sum()),in_core=int((core_depth(X)>0).sum()),secs=time.time()-t0))
    st["done"]=bool(L<1.0 or len(kinds)>=NMAX)
    json.dump(st,open(fn+".tmp","w")); os.replace(fn+".tmp",fn)
    print(f"{time.time()-t0:6.0f}s {tag} N={len(kinds)} left={L:.2f} inside={st['rounds'][-1]['inside_us']} core={st['rounds'][-1]['in_core']} m1-5={st['rounds'][-1]['spec']['rms_m1_5']:.2f} 13+={st['rounds'][-1]['spec']['rms_13up']:.2f}",flush=True)
print("DONE",tag,flush=True)

"""Fit N big masses + deep voids (core_open all-open model: anywhere 0.5..102 Gpc, inside our sphere, shell or core; no wall, no cap)
on top of a FIXED filled web. Target = R - w_web. Big bodies are coarse superstructure and are not checked against the web's spheres.
usage: python3 fit_on_web.py NPROC TAG [TAG...]   (TAG = b_random_s2, b_foam_rim_s3_core, ...) -> web_fits/<TAG>.json"""
import os,sys,json,time,numpy as np
os.environ.setdefault("OMP_NUM_THREADS","1")
from multiprocessing import Pool
from scipy.optimize import least_squares
WEBDIR=os.path.dirname(os.path.abspath(__file__)); CO="/workspace/clump_scan/core_open"; sys.path.insert(0,CO)
_c=os.getcwd(); os.chdir(CO); from co_core import *; os.chdir(_c)
AX=np.radians(15.16); TH_=np.radians(np.arange(360)); m=np.arange(1,13)
KINDS=["CCCC","CCCCVVVV","CCCCCCVVVVVV"]; NS={"CCCC":6,"CCCCVVVV":6,"CCCCCCVVVVVV":5}
SC=1/(FLOOR*np.sqrt(360))
def spec(v):
    a=2*np.mean(v[None]*np.cos(m[:,None]*(TH_-AX)),1); b=2*np.mean(v[None]*np.sin(m[:,None]*(TH_-AX)),1)
    rec=v.mean()+(a[:,None]*np.cos(m[:,None]*(TH_-AX))+b[:,None]*np.sin(m[:,None]*(TH_-AX))).sum(0)
    lo=(a[:,None]*np.cos(m[:5,None]*(TH_-AX)) if False else None)
    p15=(a[:5,None]*np.cos(m[:5,None]*(TH_-AX))+b[:5,None]*np.sin(m[:5,None]*(TH_-AX))).sum(0)
    return dict(amp=(np.hypot(a,b)/FLOOR).tolist(),sym_m=(a**2/np.maximum(a**2+b**2,1e-300)).tolist(),
                sym=float((a**2).sum()/(a**2+b**2).sum()),rms_m1_5=float(rms(p15)/FLOOR),rms_m6_12=float(rms(rec-v.mean()-p15)/FLOOR),
                rms_13up=float(rms(v-rec)/FLOOR),sym_rms=float(np.sqrt((a**2).sum()/2)/FLOOR),anti_rms=float(np.sqrt((b**2).sum()/2)/FLOOR),total=float(rms(v)/FLOOR))
def fit_one(mod,p0,T):
    lam=30.0; p=p0
    for it in range(4):
        f=lambda q: np.concatenate([(T-mod.pattern(q))*SC, lam*mod.violations(q)])
        r=least_squares(f,p,bounds=(mod.lo,mod.hi),method="trf",x_scale="jac",max_nfev=60*len(p),xtol=1e-8,ftol=1e-8); p=r.x
        if mod.violations(p).max(initial=0)<=1e-3: break
        lam*=10
    return p
def warm(kinds):
    import glob
    S=[]
    for f in (f"{CO}/raw/ao_{kinds}.json",f"/workspace/clump_scan/solution_space/all_open/raw_open/{kinds}.json"):
        if os.path.exists(f): S+=sorted([r for r in json.load(open(f)) if r.get("ok",r.get("viol",0)<=1e-3)],key=lambda r:r["left"])[:2]
    return S
def job(a):
    tag,kinds,i,T=a; mod=COAO(kinds); rng=np.random.default_rng(1000*i+len(kinds)); n=len(kinds)
    S=warm(kinds)
    if i<2 and len(S)>i: s=S[i]; p0=mod.encode(s["X"],s["M"],s["Rv"]); st="warm"
    else:
        p=np.zeros((n,4)); p[:,0]=rng.uniform(0,1,n)**2; p[:,1]=rng.uniform(-np.pi,np.pi,n); p[:,2]=np.arcsin(rng.uniform(-1,1,n)); p[:,3]=rng.uniform(0,1,n); p0=p.ravel(); st="cold"
    t0=time.time(); p=fit_one(mod,p0,T); p=fit_one(mod,p,T); X,M,Rv,Rc=mod.decode(p); wb=mod.pattern(p)
    return dict(tag=tag,kinds=kinds,start=f"{st}{i}",left=float(rms(T-wb)/FLOOR),ok=bool(mod.violations(p).max(initial=0)<=1e-3),X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),
                wb=wb.tolist(),core_depth=core_depth(X).tolist(),dist=np.linalg.norm(X,axis=1).tolist(),secs=time.time()-t0)
if __name__=="__main__":
    npr=int(sys.argv[1]); tags=sys.argv[2:]; os.makedirs(f"{WEBDIR}/web_fits",exist_ok=True); tasks=[]; WEB={}
    for tag in tags:
        fn=f"{WEBDIR}/web_fits/{tag}.json"
        if os.path.exists(fn) and json.load(open(fn)).get("done"): continue
        w=np.load(f"{WEBDIR}/final_{tag}.npz")["w"]; WEB[tag]=w; T=R-w
        tasks+=[(tag,k,i,T) for k in KINDS for i in range(NS[k])]
    tasks.sort(key=lambda t:-len(t[1])); res={t:[] for t in WEB}; t0=time.time()
    def dump(tag,done=False):
        w=WEB[tag]; T=R-w; best={}
        for r in res[tag]:
            if r["ok"] and (r["kinds"] not in best or r["left"]<best[r["kinds"]]["left"]): best[r["kinds"]]=r
        out=dict(tag=tag,done=done,web_only=spec(T),best={k:dict({q:v for q,v in b.items() if q!="wb"},spec=spec(T-np.array(b["wb"])),resid=(T-np.array(b["wb"])).tolist()) for k,b in best.items()},
                 runs=[{q:v for q,v in r.items() if q!="wb"} for r in res[tag]])
        json.dump(out,open(f"{WEBDIR}/web_fits/{tag}.json","w"))
    left_by={t:sum(1 for x in tasks if x[0]==t) for t in WEB}
    with Pool(npr) as pool:
        for r in pool.imap_unordered(job,tasks):
            res[r["tag"]].append(r); dump(r["tag"],done=len(res[r["tag"]])==left_by[r["tag"]])
            print(f"{time.time()-t0:6.0f}s {r['tag']:24s} {r['kinds']:13s} {r['start']:6s} left={r['left']:.1f} ok={r['ok']} t={r['secs']:.0f}",flush=True)
    print("done",flush=True)

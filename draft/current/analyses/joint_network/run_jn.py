"""usage: python3 run_jn.py NPROC PHASE   (PHASE A: free+cont, B: loose+web, C: strict)  -> raw/<level>_<sym>_N<N>.json"""
import os,sys,json,time,glob,numpy as np
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"): os.environ[v]="1"
from multiprocessing import Pool
from seeds import library,start
from jn_fast import JF,fit
from jn_core import *
NS=[8,12,16,20,24,32]
MODES=("warm","warm_r","cold","warm_r","cold","warm_r","cold","cold")
def kinds(N): return "C"*(N//2)+"V"*(N-N//2)
def rec(m,p,tag):
    dec=m.decode(p); X,M,Rv,Rc,dl=dec; h,a,s=m.split(dec); w=m.pattern_X(X,M); m0=JN("".join(m.k),m.level,m.mirror)
    br={}
    if m.L["cont"]: br["cont"]=float(m0.cont_pen(dec,m.L["cont"]).max(initial=0))
    if m.L["net"]: br["net"]=float(m0.net_pen(dec,m.L["net"]).max(initial=0))
    if m.L["trend"]: br["trend"]=float(m0.trend_pen(dec,m.L["trend"]).max(initial=0))
    ok=h.max(initial=0)<=1e-3 and a.max(initial=0)<=1e-3
    if not m.mirror:
        nm,npl,_,_=AC.hard_check(m.k,X,M,Rv); ok=ok and nm==0 and npl==0
    return dict(level=m.level,mirror=m.mirror,kinds="".join(m.k),start=tag,left=float(rms(R-w)/FLOOR),ok=bool(ok),
        hard=float(h.max(initial=0)),asym=float(a.max(initial=0)),soft_max=float(s.max(initial=0)),soft_rms=float(np.sqrt(np.mean(s**2))) if len(s) else 0.0,
        breakdown=br,trend_excess=m0.trend_excess(dec).tolist(),
        p=p.tolist(),X=X.tolist(),M=M.tolist(),Rv=Rv.tolist(),Rc=Rc.tolist(),delta=dl.tolist(),resid=(R-w).tolist(),amps=(amps(R-w)/FLOOR).tolist(),
        props={k:v for k,v in network_props(m.k,X,M,Rv).items()})
def job(a):
    level,mir,N,p0,tag=a; m=JF(kinds(N),level,mir); t0=time.time()
    try: p=fit(m,np.array(p0),nfev=600)
    except Exception as e:
        print("fit error",level,mir,N,tag,repr(e),flush=True); p=np.clip(np.array(p0),m.lo,m.hi)
    r=rec(m,p,tag); r["secs"]=time.time()-t0; r["start_only"]=bool(time.time()-t0<0.5); return r
def prev_best(level,sym,N,k=2):
    f=f"raw/{level}_{sym}_N{N}.json"
    if not os.path.exists(f): return []
    r=json.load(open(f)); r=[x for x in r if x["ok"]]
    # rank by total objective proxy: leftover^2 + (soft)^2 weighting by own lam
    r.sort(key=lambda x:x["left"]**2+(LEVELS[level]["lam"]*x["soft_rms"])**2*len(x["p"]))
    return [np.array(x["p"]) for x in r[:k]]
if __name__=="__main__":
    npr=int(sys.argv[1]); ph=sys.argv[2]; os.makedirs("raw",exist_ok=True); lib=library(); tasks=[]
    for si,(sym,mir) in enumerate((("m",True),("n",False))):
        for N in NS:
            rng=np.random.default_rng(1000*sum(map(ord,ph))+10*si+N)
            if ph[0]=="A":
                for lv in ("free","cont"):
                    m=JF(kinds(N),lv,mir)
                    for i,mode in enumerate(MODES): tasks.append((lv,mir,N,start(m,rng,lib,mode,N//2,N-N//2),f"{mode}{i}"))
            elif ph[0]=="B":
                for lv in ("loose","web"):
                    m=JF(kinds(N),lv,mir); P=prev_best("cont",sym,N,3)
                    for i,p in enumerate(P): tasks.append((lv,mir,N,p,f"from_cont{i}"))
                    for i,mode in enumerate(MODES[2:2+8-len(P)]): tasks.append((lv,mir,N,start(m,rng,lib,mode,N//2,N-N//2),f"{mode}{i}"))
            elif ph[0]=="C":
                m=JF(kinds(N),"strict",mir); P=prev_best("web",sym,N,3)+prev_best("loose",sym,N,2)
                for i,p in enumerate(P): tasks.append(("strict",mir,N,p,f"from_prev{i}"))
                for i in range(8-len(P)): tasks.append(("strict",mir,N,start(m,rng,lib,"cold",N//2,N-N//2),f"cold{i}"))
    tasks=[(a,b,c,d.tolist(),e) for a,b,c,d,e in tasks]
    done=set()
    for f in glob.glob("raw/*.json"):
        for x in json.load(open(f)): done.add((x["level"],x["mirror"],len(x["kinds"]),x["start"]))
    tasks=[t for t in tasks if (t[0],t[1],t[2],t[4]) not in done]
    tasks.sort(key=lambda t:-t[2]); t0=time.time(); print("tasks",len(tasks),flush=True)
    with Pool(npr,maxtasksperchild=4) as pool:
        for r in pool.imap_unordered(job,tasks):
            sym="m" if r["mirror"] else "n"; N=len(r["kinds"]); f=f"raw/{r['level']}_{sym}_N{N}.json"
            old=json.load(open(f)) if os.path.exists(f) else []; old.append(r); json.dump(old,open(f,"w"))
            print(f"{time.time()-t0:6.0f}s {r['level']:6s} {sym} N={N:2d} {r['start']:10s} left={r['left']:8.3f} ok={r['ok']} soft={r['soft_max']:.3f} {r['breakdown']} in={r['props']['n_inside']} core={r['props']['n_core']} t={r['secs']:.0f}",flush=True)
    print("done",time.time()-t0,flush=True)

"""Collect core-open results (read-only from ../../../core_open/raw) for the brief: best-vs-N curves and example layouts
(first N reaching 1x, as for the other examples), in the same example format as collected.pkl. -> core_open.pkl"""
import os,sys,json,glob,pickle,numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); CO="/workspace/clump_scan/core_open"
os.chdir(CO); sys.path.insert(0,CO)
from co_core import R,FLOOR,TH,rms,amps
runs={}
for f in glob.glob("raw/*.json"):
    for r in json.load(open(f)): runs.setdefault(r["scen"],[]).append(r)
def nb(r): return len(r["Xp"])+len(r["Xo"]) if r["scen"]=="ep" else len(r["name"])
def mode_info(v):
    a=amps(v)/FLOOR; rec=np.full(360,v.mean())
    for m in range(1,13): rec+=2*np.mean(v*np.cos(m*TH))*np.cos(m*TH)+2*np.mean(v*np.sin(m*TH))*np.sin(m*TH)
    return a,float(np.sqrt(2)*rms(v-rec)/FLOOR)
SKIP={("out",9),("ao",7)}   # undersampled classes left off the chart (as in core_open/summary.md)
out={}
for s in ("out","asym","ep","ao"):
    B={}
    for r in runs[s]:
        if not r["ok"]: continue
        n=nb(r)
        if n not in B or r["left"]<B[n]["left"]: B[n]=r
    curve={n:B[n]["left"] for n in sorted(B) if (s,n) not in SKIP}
    n1=min([n for n,v in curve.items() if v<=1.0])
    r=B[n1]; v=np.array(r["resid"]); a,hi=mode_info(v); left=float(rms(v)/FLOOR); assert abs(left-r["left"])<1e-6*max(1,left)
    if s=="ep":
        kinds=r["name"].split("+")[1]; e=dict(N=n1,kinds=kinds,X=np.array(r["Xo"]).reshape(-1,3),M=np.array(r["Mo"]),Rv=np.array(r["Rvo"]),
            partners=dict(X=np.array(r["Xp"]).reshape(-1,3),R=np.array(r["Rp"]),M=np.array(r["Mp"]),delta_us=r["delta_us"]),zoom=62)
    else: e=dict(N=n1,kinds=r["name"],X=np.array(r["X"]),M=np.array(r["M"]),Rv=np.array(r["Rv"]),zoom=45)
    e.update(resid=v,amps=a,hi=hi,left=left,n_core=int((np.array(r["core_depth"])>0).sum()),core_depth=list(np.round(np.array(r["core_depth"]),1)),dist=list(np.round(np.array(r["dist"]),1)))
    out[s]=dict(curve=curve,ex=e)
    print(s,"curve",{k:round(v,3) for k,v in curve.items()},"| example N",n1,round(left,3),"core bodies",e["n_core"],[d for d in e["core_depth"] if d>0])
pickle.dump(out,open(os.path.join(HERE,"core_open.pkl"),"wb"))

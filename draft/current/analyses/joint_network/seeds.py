"""Seed library: feasible layouts from earlier runs (read-only), plus web-like random starts."""
import os,glob,json,numpy as np
from jn_core import *
SRC=[f"{SS}/all_open/raw_open/*.json",f"{SS}/all_open/raw_smooth/*.json",f"{SS}/no_core_wall/raw_open/*.json",
     f"{SS}/no_core_wall/raw_smooth/*.json",f"{SS}/raw/*.json",f"{CS}/core_open/raw/out_*.json",f"{CS}/core_open/raw/asym_*.json",
     f"{SS}/asym/raw/*.json",f"{SS}/no_core_wall/raw_smooth_asym/*.json"]
def library():
    out=[]
    for pat in SRC:
        for f in glob.glob(pat):
            try: r=json.load(open(f))
            except Exception: continue
            r=r["runs"] if isinstance(r,dict) else r
            kinds=None
            if isinstance(r,list):
                for x in r:
                    k=x.get("kinds") or x.get("name") or os.path.basename(f)[:-5].split("_")[-1]
                    if set(k)-set("CV") or "X" not in x or x.get("viol",0)>1e-3: continue
                    out.append(dict(k=k,X=x["X"],M=x["M"],Rv=x["Rv"],left=x["left"],src=os.path.relpath(f,CS)))
    ins=f"{CS}/inside_visible/results_cap2e52.json"
    if os.path.exists(ins):
        for c in json.load(open(ins)).values():
            x=c["best_run"]; k=x["kinds"]
            if set(k)<=set("CV"): out.append(dict(k=k,X=x["X"],M=x["M"],Rv=x["Rv"],left=x["left"],src="inside_visible"))
    return out
def _order(k,X,M,Rv,nC,nV):
    """returns body lists ordered C...C V...V"""
    k=np.array(list(k)); X=np.array(X); M=np.abs(np.array(M)); Rv=np.array(Rv)
    iC=np.where(k=="C")[0][:nC]; iV=np.where(k=="V")[0][:nV]
    return [X[i] for i in iC],[M[i] for i in iC],[X[i] for i in iV],[Rv[i] for i in iV],[1.0]*len(iV)
def grow(rng,XC,MC,XV,RVv,QV,nC,nV,small=True):
    """add web-like bodies: voids next to existing voids (or in the 8-20 Gpc band), masses at void contacts"""
    XC=list(XC);MC=list(MC);XV=list(XV);RVv=list(RVv);QV=list(QV)
    while len(XV)<nV:
        r=rng.uniform(0.5,1.0) if small else rng.uniform(1.5,3.5)
        if XV and rng.uniform()<0.7 and not small:
            j=rng.integers(len(XV)); v=rng.normal(size=3); v/=np.linalg.norm(v); x=np.array(XV[j])+v*(RVv[j]+r)*1.05
        else:
            v=rng.normal(size=3); v/=np.linalg.norm(v); x=v*(rng.uniform(20,30) if small else rng.uniform(8,20))
        d=np.linalg.norm(x)
        if d<r+0.6 or abs(d-RV)<r+0.15: continue
        XV.append(x); RVv.append(r); QV.append(0.02 if small else rng.uniform(0.3,1))
    while len(XC)<nC:
        if small:
            v=rng.normal(size=3); v/=np.linalg.norm(v); x=v*rng.uniform(20,30)
        elif len(XV)>=2:
            a,b=rng.choice(len(XV),2,replace=False); A=np.array(XV[a]); B=np.array(XV[b])
            x=A+(B-A)*RVv[a]/(RVv[a]+RVv[b])+rng.normal(size=3)*0.5
        else:
            v=rng.normal(size=3); v/=np.linalg.norm(v); x=v*rng.uniform(8,20)
        d=np.linalg.norm(x)
        if d<1.5 or abs(d-RV)<1.0: continue
        XC.append(x); MC.append(10**(49 if small else rng.uniform(50.5,52)))
    return XC,MC,XV,RVv,QV
def start(m,rng,lib,mode,nC,nV):
    if mode=="cold":
        XC,MC,XV,RVv,QV=grow(rng,[],[],[],[],[],nC,nV,small=False)
    else:
        cand=[s for s in lib if s["k"].count("C")<=nC and s["k"].count("V")<=nV and s["k"].count("C")>=1]
        cand.sort(key=lambda s:s["left"])
        top=cand[:8]
        s=top[rng.integers(len(top))] if (mode=="warm_r" and top) else (top[0] if top else None)
        if s is None: return start(m,rng,lib,"cold",nC,nV)
        XC,MC,XV,RVv,QV=_order(s["k"],s["X"],s["M"],s["Rv"],nC,nV)
        XC,MC,XV,RVv,QV=grow(rng,XC,MC,XV,RVv,QV,nC,nV,small=True)
    X=np.array(XC+XV); M=np.array(MC+[1.0]*nV); Rv=np.array([1.0]*nC+RVv); q=np.array([1.0]*nC+QV)
    p=m.encode(X,M,Rv,q)
    if mode=="warm_r": p=p+rng.normal(0,0.02,len(p))*(m.hi-m.lo)*0.05
    return np.clip(p,m.lo,m.hi)

"""Cluster multi-start solutions into distinct basins, summarise, make charts. -> results.json, summary tables, PNGs."""
import json, glob, os, numpy as np
from scipy.optimize import linear_sum_assignment
import sys; sys.path.insert(0,"."); from ss_core import *
TH_L=[89,50,27,10]; ACTIVE=1.0      # object counts as active if its own pattern is >= 1x floor (rms)
POS_TOL=2.0; LOGM_TOL=0.25; RV_TOL=0.5
def mirror(X):
    a=np.radians(2*AXIS); c,s=np.cos(a),np.sin(a); x,y=X[:,0],X[:,1]
    return np.c_[c*x+s*y,s*x-c*y,X[:,2]]
def canon(run,kinds):
    X=np.array(run["X"]); X[:,2]=np.abs(X[:,2]); M=np.array(run["M"]); Rv=np.array(run["Rv"]); k=np.array(list(kinds))
    act=np.zeros(len(k),bool)
    for i in range(len(k)):
        w=M[i]*unit_field(X[i:i+1])[0] if k[i]!="S" else shallow_field(X[i:i+1],Rv[i])[0]
        act[i]=rms(w)/FLOOR>=ACTIVE
    size=np.where(k=="C",np.log10(np.abs(M)+1),Rv)
    return dict(X=X[act],k=k[act],size=size[act],left=run["left"])
def same(a,b):
    if sorted(a["k"])!=sorted(b["k"]): return False
    for Xb in (b["X"],mirror(b["X"])):
        ok=True
        for t in set(a["k"]):
            ia=np.where(a["k"]==t)[0]; ib=np.where(b["k"]==t)[0]
            D=np.linalg.norm(a["X"][ia][:,None]-Xb[ib][None],axis=2); r,c=linear_sum_assignment(D)
            ds=np.abs(a["size"][ia][r]-b["size"][ib][c]).max(); tol=LOGM_TOL if t=="C" else RV_TOL
            if D[r,c].max()>POS_TOL or ds>tol: ok=False; break
        if ok: return True
    return False
def mix_of(k):
    s=set(k)
    return {"C":"masses only","V":"deep voids only","S":"shallow voids only"}[k[0]] if len(s)==1 else ("masses + deep voids" if "V" in s else "masses + shallow voids")
res={}
for f in sorted(glob.glob("raw/*.json")):
    d=json.load(open(f)); k=d["kinds"]; runs=[r for r in d["runs"] if r["viol"]<=1e-3]; L=np.array(sorted(r["left"] for r in runs))
    runs40=[r for r in d["runs"][:40] if r["viol"]<=1e-3]   # common 40-start budget for fair distinct counts
    good=sorted([r for r in runs40 if r["left"]<TH_L[0]],key=lambda r:r["left"]); reps=[]
    for r in good:
        c=canon(r,k)
        for rep in reps:
            if same(rep,c): rep["n"]+=1; break
        else: c["n"]=1; reps.append(c)
    res[k]=dict(kinds=k,N=len(k),n_C=k.count("C"),n_V=k.count("V"),n_S=k.count("S"),mix=mix_of(k),params=4*len(k),
        starts=d["nstarts"],feasible=len(runs),best=float(L[0]),p10=float(np.percentile(L,10)),median=float(np.median(L)),
        frac_below={str(t):float((L<t).mean()) for t in TH_L},n_runs_for_distinct=len(runs40),
        distinct_below={str(t):int(sum(rep["left"]<t for rep in reps)) for t in TH_L},
        basins=[dict(left=rep["left"],hits=rep["n"],n_active=len(rep["k"]),kinds="".join(rep["k"]),X=rep["X"].tolist(),size=rep["size"].tolist()) for rep in reps],
        best_amps=min(runs,key=lambda r:r["left"])["amps"],best_run=min(runs,key=lambda r:r["left"]))
json.dump(res,open("results.json","w"))

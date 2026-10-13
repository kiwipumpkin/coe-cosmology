"""Distinct-basin counts, caps, multipoles for inside-sphere fits -> results_TAG.json"""
import json, glob, os, sys, numpy as np
from scipy.optimize import linear_sum_assignment
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); from inside_core import *
from ss_core import AXIS
TH_L=[50,27,10,1]; ACTIVE=1.0; POS_TOL=2.0; LOGM_TOL=0.25; RV_TOL=0.5
def mirror(X):
    a=np.radians(2*AXIS); c,s=np.cos(a),np.sin(a); x,y=X[:,0],X[:,1]; return np.c_[c*x+s*y,s*x-c*y,X[:,2]]
def canon(run,kinds):
    X=np.array(run["X"]); X[:,2]=np.abs(X[:,2]); M=np.array(run["M"]); Rv=np.array(run["Rv"]); k=np.array(list(kinds))
    act=np.array([rms(M[i]*unit_field(X[i:i+1])[0])/FLOOR>=ACTIVE for i in range(len(k))])
    size=np.where(k=="C",np.log10(np.abs(M)+1),Rv); return dict(X=X[act],k=k[act],size=size[act],left=run["left"])
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
    s=set(k); return {"C":"masses only","V":"deep voids only"}[k[0]] if len(s)==1 else "masses + deep voids"
def analyze(tag):
    res={}
    for f in glob.glob(f"raw_{tag}/*.json"):
        d=json.load(open(f)); k=d["kinds"]; mhi=d["mhi"]; runs=[r for r in d["runs"] if r["viol"]<=1e-3]; L=np.array(sorted(r["left"] for r in runs))
        good=sorted([r for r in runs if r["left"]<TH_L[0]],key=lambda r:r["left"]); reps=[]
        for r in good:
            c=canon(r,k)
            for rep in reps:
                if same(rep,c): rep["n"]+=1; break
            else: c["n"]=1; reps.append(c)
        b=min(runs,key=lambda r:r["left"]); kk=np.array(list(k)); M=np.array(b["M"]); Rv=np.array(b["Rv"]); Rc=np.array(b["Rc"]); X=np.array(b["X"]); dd=np.linalg.norm(X,axis=1)
        res[k]=dict(kinds=k,N=len(k),n_C=k.count("C"),n_V=k.count("V"),mix=mix_of(k),starts=d["nstarts"],feasible=len(runs),best=float(L[0]),median=float(np.median(L)),
            distinct_below={str(t):int(sum(rep["left"]<t for rep in reps)) for t in TH_L},
            best_caps=dict(mass_at_cap=int(((kk=="C")&(M>=0.98*mhi)).sum()),void_at_4Gpc=int(((kk=="V")&(Rv>=3.98)).sum()),
                           at_ring_margin=int((dd+Rc>=RV-0.1-0.05).sum()),at_us=int((dd-Rc<=0.1+0.05).sum())),
            best_masses_kg=float(M[kk=="C"].sum()),best_void_deficit_kg=float(M[kk=="V"].sum()),best_net_kg=float(M.sum()),
            best_dist_gpc=dd.tolist(),best_amps=b["amps"],best_left_surroundcomp=b["left_surroundcomp"],best_left_nocomp=b["left_nocomp"],best_run=b)
    json.dump(res,open(f"results_{tag}.json","w")); return res
if __name__=="__main__":
    for t in sys.argv[1:]:
        r=analyze(t); print(t,len(r),"classes")

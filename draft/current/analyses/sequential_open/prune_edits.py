"""Prune a minimal-edit end state: undo edits one at a time (restore the body as grown) whenever the leftover stays <= LIM (default 1x)
and the restored body overlaps nothing. Several passes; order = smallest pattern change first. Writes minedit2/<tag>_pruned.json (same schema)."""
import json,sys,numpy as np
from scipy.spatial import cKDTree
from common import *
from fit_on_web import spec
K=seq_core.unit_pattern_fast; GAP=0.02; LIM=float(sys.argv[2]) if len(sys.argv)>2 else 1.0
js=sys.argv[1]; d=json.load(open(js)); base=d["tag"].split("_clr")[0]
f=np.load(f"{OUT}/final_{base}.npz"); X0=f["Xp"].astype(float); r0=f["rp"].astype(float); M0=f["Mp"].astype(float); isV=f["kp"]=="void"; n=len(X0)
idx=np.array(d["idx"],int); X=X0.copy(); r=r0.copy(); M=M0.copy(); X[idx]=d["X_new"]; r[idx]=d["r_new"]; M[idx]=d["M_new"]
W=np.zeros(360)
for i in range(0,n,3000): W+=M[i:i+3000]@K(X[i:i+3000])
L0=rms(R-W)/FLOOR; print(js,"start",round(L0,3),flush=True)
lim=max(LIM,L0)
Pn=K(X[idx]); Po=K(X0[idx]); dW=M0[idx][:,None]*Po-M[idx][:,None]*Pn     # change of pattern if body k is restored
undone=np.zeros(len(idx),bool)
for ps in range(4):
    order=np.argsort(np.linalg.norm(dW,axis=1)); nu=0
    alive=r>0; tr=cKDTree(X[alive]); ra=r[alive]; ai=np.where(alive)[0]
    for k in order:
        if undone[k]: continue
        Wt=W+dW[k]
        if rms(R-Wt)/FLOOR>lim: continue
        i=idx[k]; nb=tr.query_ball_point(X0[i],r0[i]+ra.max()+GAP)
        ok=all((ai[j]==i) or (np.linalg.norm(X[ai[j]]-X0[i])>=r[ai[j]]+r0[i]+GAP) for j in nb)
        if not ok: continue
        W=Wt; X[i]=X0[i]; r[i]=r0[i]; M[i]=M0[i]; undone[k]=True; nu+=1
        alive=r>0; tr=cKDTree(X[alive]); ra=r[alive]; ai=np.where(alive)[0]
    print(" pass",ps,"undone",nu,"left",round(rms(R-W)/FLOOR,3),flush=True)
    if nu==0: break
keep=~undone; e=idx[keep]
for key in ("is_void","near_ring","deleted","moved","resized","a_new","X_old","X_new","r_old","r_new","M_old","M_new","move","dist_us_old","dist_us_new","dist_ring_old","dist_ring_new"):
    if key in d: d[key]=list(np.array(d[key],dtype=object)[keep]) if not isinstance(d[key][0],(int,float,bool)) else np.array(d[key])[keep].tolist()
for key in ("X_old","X_new"): d[key]=[list(map(float,x)) for x in d[key]]
d["idx"]=e.tolist(); res=R-W; d["resid"]=res.tolist(); d["after"]=spec(res); d["left"]=float(rms(res)/FLOOR); d["reached"]=bool(d["left"]<=1)
de=np.array(d["deleted"],bool); d["n_deleted"]=int(de.sum()); d["n_moved"]=int(np.array(d["moved"],bool).sum()); d["n_resized"]=int(np.array(d["resized"],bool).sum())
d["n_edited"]=int(len(e)); d["frac_edited"]=float(len(e)/n); d["n_near_ring_edited"]=int(np.array(d["near_ring"],bool).sum()); d["tag"]=d["tag"]+"_pruned"; d["pruned_undone"]=int(undone.sum())
d["net_mass_change"]=float(M[e].sum()-M0[e].sum())
d["curve"]=d["curve"]+[(int(len(e)),float(d["left"]))]
json.dump(d,open(js.replace(".json","_pruned.json"),"w")); print("DONE",d["tag"],"edited",len(e),"deleted",d["n_deleted"],"left",round(d["left"],3),flush=True)

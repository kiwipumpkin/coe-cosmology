"""Shared features of good solutions (distinct basins below 50x and 27x)."""
import json, numpy as np, sys; sys.path.insert(0,"."); from ss_core import *
res=json.load(open("results.json")); out={}
def feat(thr,filt):
    rows=[]
    for k,r in res.items():
        if not filt(r): continue
        for b in r["basins"]:
            if b["left"]>=thr: continue
            X=np.array(b["X"]).reshape(-1,3); kk=np.array(list(b["kinds"])); s=np.array(b["size"])
            d=np.linalg.norm(X,axis=1); l=np.degrees(np.arctan2(X[:,1],X[:,0]))%360; lat=np.degrees(np.arcsin(X[:,2]/d))
            for i in range(len(kk)): rows.append(dict(cls=k,kind=kk[i],d=d[i],l=l[i],lat=lat[i],size=s[i],dcoe=np.linalg.norm(X[i]-C_G)))
    if not rows: return None
    f={}
    for t in "CVS":
        R_=[x for x in rows if x["kind"]==t]
        if not R_: continue
        d=np.array([x["d"] for x in R_]); lat=np.array([x["lat"] for x in R_]); l=np.array([x["l"] for x in R_]); s=np.array([x["size"] for x in R_])
        hist=np.histogram(l,bins=np.arange(0,361,30))[0]
        f[t]=dict(n=len(R_),dist_median=float(np.median(d)),dist_p10_p90=[float(np.percentile(d,10)),float(np.percentile(d,90))],
                  frac_within_20gpc=float((d<20).mean()),frac_lat_below15=float((lat<15).mean()),lat_median=float(np.median(lat)),
                  l_hist_30deg=hist.tolist(),frac_l_270_300=float(((l>=270)&(l<300)).mean()),
                  frac_at_upper_size_bound=float((s>=(np.log10(2e52)-0.01 if t=="C" else 3.99)).mean()),
                  size_median=float(np.median(s)))
    return f
for thr in (50,27):
    out[f"below{thr}_all_with_masses"]=feat(thr,lambda r:r["n_C"]>0)
    out[f"below{thr}_masses_only"]=feat(thr,lambda r:r["mix"]=="masses only")
    out[f"below{thr}_masses_deep"]=feat(thr,lambda r:r["mix"]=="masses + deep voids")
# total masses per good basin from best runs
tot=[]
for k,r in res.items():
    b=r["best_run"]; M=np.array(b["M"]); kk=np.array(list(k))
    tot.append(dict(cls=k,left=r["best"],M_clusters=float(M[kk=="C"].sum()),M_deepvoids=float(M[kk=="V"].sum()),
        max_cluster=float(M[kk=="C"].max()) if (kk=="C").any() else 0.0))
out["best_run_totals"]=sorted(tot,key=lambda x:x["left"])
json.dump(out,open("features.json","w"),indent=1)

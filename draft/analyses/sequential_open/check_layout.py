"""Legality check of a minimal-edit end state: overlaps among surviving bodies, moved bodies within 3 Gpc of the ring, region."""
import numpy as np,json,sys
from scipy.spatial import cKDTree
from common import *
def check(js,tag):
    d=json.load(open(js)); f=np.load(f"final_{tag}.npz"); X=f["Xp"].astype(float).copy(); r=f["rp"].astype(float).copy()
    i=np.array(d["idx"],int); X[i]=np.array(d["X_new"]); r[i]=np.array(d["r_new"]); alive=r>0
    Xa,ra=X[alive],r[alive]; t=cKDTree(Xa); p=np.array(list(t.query_pairs(2*ra.max())))
    pen=ra[p[:,0]]+ra[p[:,1]]-np.linalg.norm(Xa[p[:,0]]-Xa[p[:,1]],axis=1); ov=pen>1e-3
    edited=np.zeros(len(X),bool); edited[i]=True; ea=edited[alive]
    ov_e=ov&(ea[p[:,0]]|ea[p[:,1]])
    mv=np.array(d.get("moved",[True]*len(i)),bool); dr=np.sqrt((np.hypot(X[i,0],X[i,1])-RMIN)**2+X[i,2]**2)-r[i]
    ringbad=int(((dr<3-1e-3)&mv&(r[i]>0)).sum()); s=np.linalg.norm(X[i],axis=1); reg=int((((s+r[i])>30+1e-3)|((s-r[i])<0.5-1e-3))[(r[i]>0)&mv].sum())
    return dict(overlap_pairs=int(ov.sum()),overlap_pairs_involving_edited=int(ov_e.sum()),max_overlap_gpc=float(pen[ov].max()) if ov.any() else 0.0,moved_within_3gpc=ringbad,outside_region=reg)
if __name__=="__main__":
    for js in sys.argv[1:]:
        tag=json.load(open(js))["tag"].split("_clr")[0].replace("_nr","")
        print(js,check(js,tag))

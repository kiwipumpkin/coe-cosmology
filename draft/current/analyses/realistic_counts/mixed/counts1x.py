"""At 1x: smallest-change box solution (largest L2 damping that still reaches 1x) for given X; composition per cell and per lobe."""
import os, sys, json, numpy as np
from mix_core import *
from scan import CASES
RHO_G=RHO*GPC**3                      # kg per Gpc^3 at mean density
E_c=DV_FRAC*RHO_G*S_GPC**3/NC_PER_CELL; D_v=DV_FRAC*RHO_G*S_GPC**3
def analyse(gs,ds):
    d=np.concatenate(ds); vol=np.concatenate([np.full(len(g.P),g.H**3) for g in gs]); P=np.concatenate([g.P for g in gs])
    a,b,q=composition(d); conv=np.abs(a)*n_c0*vol                      # clusters converted (= voids converted)
    act=np.abs(d)>0.005
    return dict(M_abs=float((np.abs(d)*rho_gpc(P)*(GPC**3)*vol).sum()),converted_clusters=float(conv.sum()),converted_voids=float(conv.sum()),
        cluster_only_count=float((np.abs(d)*RHO_G*vol).sum()/E_c),void_only_count=float((np.abs(d)*RHO_G*vol).sum()/D_v),
        objects_present_active=float(((n_c0+n_v0)*vol[act]).sum()),vol_active=float(vol[act].sum()),peak=float(np.abs(d).max()),
        q_min=float(q.min()),q_max=float(q.max()),voidfill_min=float((FV0*(1+b)).min()),voidfill_max=float((FV0*(1+b)).max()))
if __name__=="__main__":
    case=sys.argv[1]; Xs=[float(x) for x in sys.argv[2].split(",")]; gs=CASES[case](); out={}
    for X in Xs:
        for lam in [1000,100,10,1,0]:
            ds,L=solve_box(gs,X,1.0,lam)
            if L<=1.0: break
        r=analyse(gs,ds); r.update(X=X,lam=lam,left=L); out[str(X)]=r
        np.savez(f"best1x_{case}_X{X}.npz",*ds); print(case,json.dumps(r),flush=True)
        json.dump(out,open(f"counts1x_{case}.json","w"),indent=1)

import os,json,numpy as np
from scipy.optimize import lsq_linear
from cb_core import *
sys.path.insert(0,"/workspace/clump_scan/solution_space"); _c=os.getcwd(); import ss_core; os.chdir(_c)
A=0.10; OFF=1.5; MMAX=2e52; NB=120
def cand(rng,n=600):
    out=[]
    while len(out)<n:
        x=rng.uniform(-30,30,3); s=np.linalg.norm(x)
        if s>30 or s<1: continue
        if np.linalg.norm(x-C_G)<A_IN_G+0.5: continue
        if np.hypot(np.hypot(x[0],x[1])-RV,x[2])<OFF: continue   # >= 1.5 Gpc from every ring point
        out.append(x)
    return np.array(out)
def fit(target,rng):
    X=cand(rng); U=ss_core.unit_field(X); U-=U.mean(1,keepdims=True)
    WL=8.0   # weight on m<=12 vs the rest: minimise the full leftover, with m<=12 pushed below 1x
    lo=lambda v: lowpass(v); L=np.array([WL*lo(u)+(u-lo(u)) for u in U]).T; t=-(WL*lo(target)+(target-lo(target)))
    sc=MMAX
    r=lsq_linear(L*sc,t,bounds=(-1,1),lsmr_tol='auto',max_iter=5000)
    keep=np.argsort(-np.abs(r.x))[:NB]
    r2=lsq_linear(L[:,keep]*sc,t,bounds=(-1,1),max_iter=5000)
    Mk=r2.x*sc; w=Mk@U[keep]
    return X[keep],Mk,w
pts=points()
rng=np.random.default_rng(7); out={}; table=[]
for mdl in ["chess","einasto"]:
    for i in range(9):
        Q=Q_SG if i==0 else rand_Q(rng); ph=(0,0,0) if i==0 else rng.uniform(0,1,3)
        w,wus=pattern(mdl,A,Q,ph,pts=pts)
        Xb,Mb,wb=fit(R+w,rng); tot=R+w+wb
        Xa,Ma,wa=fit(w,rng); tot_a=w+wa        # checkerboard alone, cancelled
        row=dict(model=mdl,run=i,orient="supergalactic (Einasto)" if i==0 else "random",
            own=rms(w)/FLOOR,own_m1_12=band(w,1,12),own_13=band(w,13,999),own_sym=sym_share(w),own_dipole_us=rms(wus-wus.mean())/FLOOR,
            comb_before=rms(R+w)/FLOOR,comb_before_m1_12=band(R+w,1,12),comb_before_13=band(R+w,13,999),comb_before_sym=sym_share(R+w),
            after=rms(tot)/FLOOR,after_m1_12=band(tot,1,12),after_m1_5=band(tot,1,5),after_m6_12=band(tot,6,12),after_13=band(tot,13,999),after_sym=sym_share(tot),
            n_bodies=int((np.abs(Mb)>1e47).sum()),n_pos=int((Mb>1e47).sum()),n_neg=int((Mb<-1e47).sum()),Mabs_max=float(np.abs(Mb).max()),
            min_ring_dist=float(np.min(np.hypot(np.hypot(Xb[:,0],Xb[:,1])-RV,Xb[:,2]))),
            alone_after_m1_12=band(tot_a,1,12),alone_after_13=band(tot_a,13,999),alone_after_sym=sym_share(tot_a))
        print({k:(round(v,3) if isinstance(v,float) else v) for k,v in row.items()},flush=True); table.append(row)
        if i==0: out[mdl]=dict(w=w,wus=wus,Xb=Xb,Mb=Mb,wb=wb,tot=tot,wa=wa,tot_a=tot_a,Q=Q)
np.savez("main_runs.npz",**{f"{m}_{k}":v for m in out for k,v in out[m].items()},R=R,FLOOR=FLOOR)
json.dump(table,open("runs.json","w"),indent=1)

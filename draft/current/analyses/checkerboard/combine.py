import json,numpy as np
exec(open("main_fit.py").read().split("pts=points()")[0])
d=dict(np.load("main_fit.npz"))
best=None
for sd in [11,31,41,51]:
    Xr_,Mr_,wr_,_=fit(R.copy(),np.random.default_rng(sd)); sc=band(R+wr_,13,999)
    print("seed",sd,band(R+wr_,1,12),sc,flush=True)
    if band(R+wr_,1,12)<0.9 and (best is None or sc<best[0]): best=(sc,Xr_,Mr_,wr_)
_,Xr,Mr,wr=best                    # set 1: bodies that cancel the SR+central-mass residual
res={"none":dict(after_m1_12=band(R+wr,1,12),after_13=band(R+wr,13,999),bodies_13=band(wr,13,999),n=int((abs(Mr)>1e47).sum()))}
out=dict(R=R,FLOOR=FLOOR,Xr=Xr,Mr=Mr,wr=wr)
for mdl in ["chess","einasto"]:
    w=d[f"{mdl}_w"]; wa=d[f"{mdl}_wa"]; Xa=d[f"{mdl}_Xa"]; Ma=d[f"{mdl}_Ma"]
    tot=R+w+wr+wa; Xall=np.r_[Xr,Xa]; Mall=np.r_[Mr,Ma]
    ring=np.hypot(np.hypot(Xall[:,0],Xall[:,1])-RV,Xall[:,2])
    res[mdl]=dict(A=0.10,own=rms(w)/FLOOR,own_m1_12=band(w,1,12),own_m1_5=band(w,1,5),own_13=band(w,13,999),own_m1_only=band(w,1,1),own_sym=sym_share(w),
        before=rms(R+w)/FLOOR,before_m1_12=band(R+w,1,12),before_13=band(R+w,13,999),before_sym=sym_share(R+w),
        after=rms(tot)/FLOOR,after_m1_5=band(tot,1,5),after_m6_12=band(tot,6,12),after_m1_12=band(tot,1,12),after_13=band(tot,13,999),
        after_sym=sym_share(tot),after_sym_13=sym_share(tot,13),bodies_13=band(wr+wa,13,999),
        n_bodies=int((abs(Mall)>1e47).sum()),n_pos=int((Mall>1e47).sum()),n_neg=int((Mall<-1e47).sum()),
        min_ring_dist_Gpc=float(ring[abs(Mall)>1e47].min()),Mabs_max=float(abs(Mall).max()),
        n_inside_sphere=int(((np.linalg.norm(Xall,axis=1)<RV)&(abs(Mall)>1e47)).sum()))
    out.update({f"{mdl}_w":w,f"{mdl}_wa":wa,f"{mdl}_tot":tot,f"{mdl}_X":Xall,f"{mdl}_M":Mall,f"{mdl}_wus":d[f"{mdl}_wus"]})
    print(mdl,{k:round(float(v),3) if isinstance(v,(float,np.floating)) else v for k,v in res[mdl].items()})
print(res["none"])
np.savez("final.npz",**out); json.dump(res,open("final.json","w"),indent=1,default=float)

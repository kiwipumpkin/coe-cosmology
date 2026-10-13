import os,json,numpy as np
from scipy.optimize import lsq_linear
from cb_core import *
sys.path.insert(0,"/workspace/clump_scan/solution_space"); _c=os.getcwd(); import ss_core; os.chdir(_c)
A=0.10; OFF=1.5; MMAX=2e52
def cand(rng,n):
    out=[]
    while len(out)<n:
        x=rng.uniform(-30,30,3); s=np.linalg.norm(x)
        if s>30 or s<1 or np.linalg.norm(x-C_G)<A_IN_G+0.5 or np.hypot(np.hypot(x[0],x[1])-RV,x[2])<OFF: continue
        out.append(x)
    return np.array(out)
def fit(target,rng):
    for WL,NB,NC in [(8,120,600),(25,150,900),(80,200,1200),(250,250,1500)]:
        X=cand(rng,NC); U=ss_core.unit_field(X); U-=U.mean(1,keepdims=True)
        L=np.array([WL*lowpass(u)+(u-lowpass(u)) for u in U]).T; t=-(WL*lowpass(target)+(target-lowpass(target)))
        r=lsq_linear(L*MMAX,t,bounds=(-1,1),max_iter=5000); keep=np.argsort(-np.abs(r.x))[:NB]
        r2=lsq_linear(L[:,keep]*MMAX,t,bounds=(-1,1),max_iter=5000); M=r2.x*MMAX; w=M@U[keep]
        print("  WL",WL,"NB",NB,"m1-12",band(target+w,1,12),"13+",band(target+w,13,999),flush=True)
        if band(target+w,1,12)<0.9: break
    return X[keep],M,w,WL
pts=points(); rng=np.random.default_rng(11); out={}; info={}
for mdl in ["chess","einasto"]:
    w,wus=pattern(mdl,A,pts=pts)
    Xb,Mb,wb,WL=fit(R+w,rng); Xa,Ma,wa,_=fit(w,rng)
    out.update({f"{mdl}_{k}":v for k,v in dict(w=w,wus=wus,Xb=Xb,Mb=Mb,wb=wb,Xa=Xa,Ma=Ma,wa=wa).items()})
    tot=R+w+wb; info[mdl]=dict(own=rms(w)/FLOOR,own_m1_12=band(w,1,12),own_13=band(w,13,999),own_sym=sym_share(w),own_sym_13=sym_share(w,13),
       before=rms(R+w)/FLOOR,before_m1_12=band(R+w,1,12),before_13=band(R+w,13,999),before_sym=sym_share(R+w),
       after=rms(tot)/FLOOR,after_m1_5=band(tot,1,5),after_m6_12=band(tot,6,12),after_m1_12=band(tot,1,12),after_13=band(tot,13,999),after_sym=sym_share(tot),after_sym_13=sym_share(tot,13),
       n_bodies=int((np.abs(Mb)>1e47).sum()),n_pos=int((Mb>1e47).sum()),n_neg=int((Mb<-1e47).sum()),Mabs_max=float(np.abs(Mb).max()),Mabs_median=float(np.median(np.abs(Mb[np.abs(Mb)>1e47]))),
       min_ring_dist=float(np.min(np.hypot(np.hypot(Xb[:,0],Xb[:,1])-RV,Xb[:,2]))),WL=WL,
       alone_after_m1_12=band(w+wa,1,12),alone_after_13=band(w+wa,13,999),alone_after_sym=sym_share(w+wa))
    print(mdl,info[mdl],flush=True)
# R alone fitted with the same bodies machinery (no background) for reference
Xr,Mr,wr,_=fit(R.copy(),rng); out.update(R_wb=wr); info["none"]=dict(after_m1_12=band(R+wr,1,12),after_13=band(R+wr,13,999),after=rms(R+wr)/FLOOR,after_sym=sym_share(R+wr))
print("none",info["none"])
np.savez("main_fit.npz",R=R,FLOOR=FLOOR,**out); json.dump(info,open("main_fit.json","w"),indent=1)

import json,numpy as np
import sys; sys.path.insert(0,"/workspace/clump_scan/checkerboard")
exec(open("/workspace/clump_scan/checkerboard/main_fit.py").read().split("pts=points()")[0])
from ri_core import delta_r, pattern_r
cb=np.load("/workspace/clump_scan/checkerboard/final.npz"); Xr,Mr,wr=cb["Xr"],cb["Mr"],cb["wr"]   # 150 bodies cancelling the 237x residual (same set as checkerboard)
pts=points(1.0,LAM/6,sphere_only=True); rng=np.random.default_rng(5); res={}; out=dict(R=R,FLOOR=FLOOR,Xr=Xr,Mr=Mr,wr=wr)
for name,lam in [("einasto171",LAM),("bao150",0.150)]:
    w,wus=pattern_r(0.10,lam,pts=pts)
    ens=[]
    for phi in np.linspace(0,2*np.pi,8,endpoint=False)[1:]:
        we,_=pattern_r(0.10,lam,phi,pts=pts); ens.append([rms(we)/FLOOR,band(we,1,12),band(we,13,999),sym_share(we),sym_share(we,1,12)])
    ens=np.array(ens)
    Xa,Ma,wa,_=fit(w,rng); tot=R+w+wr+wa; Xall=np.r_[Xr,Xa]; Mall=np.r_[Mr,Ma]; s=np.abs(Mall)>1e47
    ring=np.hypot(np.hypot(Xall[:,0],Xall[:,1])-RV,Xall[:,2])
    res[name]=dict(lam_Mpc=lam*1000,A=0.10,own=rms(w)/FLOOR,own_m1_5=band(w,1,5),own_m6_12=band(w,6,12),own_m1_12=band(w,1,12),own_13=band(w,13,999),own_sym=sym_share(w),own_sym_1_12=sym_share(w,1,12),own_sym_13=sym_share(w,13),
      phase_ens=dict(own=[ens[:,0].min(),ens[:,0].max()],m1_12=[ens[:,1].min(),ens[:,1].max()],m13=[ens[:,2].min(),ens[:,2].max()],sym=[ens[:,3].min(),ens[:,3].max()],sym_1_12=[ens[:,4].min(),ens[:,4].max()]),
      before=rms(R+w)/FLOOR,before_m1_12=band(R+w,1,12),before_13=band(R+w,13,999),before_sym=sym_share(R+w),
      after=rms(tot)/FLOOR,after_m1_5=band(tot,1,5),after_m6_12=band(tot,6,12),after_m1_12=band(tot,1,12),after_13=band(tot,13,999),after_sym=sym_share(tot),after_sym_13=sym_share(tot,13),
      alone_after_13=band(w+wa,13,999),bodies_13=band(wr+wa,13,999),
      n_bodies=int(s.sum()),n_pos=int((Mall>1e47).sum()),n_neg=int((Mall<-1e47).sum()),min_ring_dist_Gpc=float(ring[s].min()),Mabs_max=float(np.abs(Mall).max()))
    out.update({f"{name}_{k}":v for k,v in dict(w=w,wa=wa,tot=tot,X=Xall,M=Mall).items()})
    print(name,json.dumps(res[name],default=float),flush=True)
np.savez("final.npz",**out); json.dump(res,open("final.json","w"),indent=1,default=float)

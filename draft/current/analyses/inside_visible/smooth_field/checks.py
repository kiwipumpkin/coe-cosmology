import json, numpy as np
from solve import *
from lcdm import T_nw, A as AN, h, ns
res={}
# 1) where the 1x field lives (H=0.3, r_max 13.8)
g=Grid(0.3)
r=np.linalg.norm(g.P,axis=1); zabs=np.abs(g.P[:,2]); lon=np.degrees(np.arctan2(g.P[:,1],g.P[:,0]))%360
for s in [0,1,2,4]:
    d=np.load(f"delta_H0.3_s{s}.npy"); w=d**2/np.sum(d**2)
    bins=[0.5,4,8,11,12.5,13.8]; prof=[float(w[(r>=a)&(r<b)].sum()) for a,b in zip(bins[:-1],bins[1:])]
    res[f"s{s}"]=dict(power_frac_by_radius=dict(zip([f"{a}-{b}" for a,b in zip(bins[:-1],bins[1:])],prof)),
        frac_within_2Gpc_of_plane=float(w[zabs<2].sum()),frac_r_gt_11=float(w[r>11].sum()),
        mass_moved_kg=float(np.abs(g.m*d).sum()/2),mass_moved_frac=float(np.abs(g.m*d).sum()/2/M_VIS))
    print(s,res[f"s{s}"],flush=True)
# 2) ring-margin sensitivity
for rmax in [13.3,12.3,10.3]:
    gg=Grid(0.3,r_max=rmax)
    for s in [0,2]:
        S=Solver(gg,s); _,sc=S.to_target(1.0); res[f"rmax{rmax}_s{s}"]=sc; print("rmax",rmax,"s",s,{k:round(v,4) for k,v in sc.items()},flush=True)
# 3) what a random LCDM field (z=0 linear, scales >= grid) does through this kernel
n=g.n; H=g.H; k1=2*np.pi*np.fft.fftfreq(n,d=H)            # 1/Gpc
KX,KY,KZ=np.meshgrid(k1,k1,np.fft.rfftfreq(n,d=H)*2*np.pi,indexing="ij"); K=np.sqrt(KX**2+KY**2+KZ**2)
kh=K/1000/h; P=np.where(K>0,AN*kh**ns*T_nw(np.maximum(kh,1e-8))**2,0)/h**3*1e-9    # P in Gpc^3 (P[(Mpc/h)^3] -> Gpc^3)
rng=np.random.default_rng(3); L=[];sig=[]
for i in range(20):
    wn=np.fft.rfftn(rng.normal(size=(n,n,n)))
    f=np.fft.irfftn(wn*np.sqrt(P/H**3),s=(n,n,n)); d=f[g.mask]
    y=g.A@d; L.append(float(np.sqrt((y**2).sum()))); sig.append(float(d.std()))
res["lcdm_random"]=dict(pattern_rms_over_floor=L,delta_rms_grid=sig,median_pattern=float(np.median(L)))
print("LCDM random fields: pattern rms/floor median",np.median(L),"range",min(L),max(L),"delta rms on 0.3Gpc grid",np.median(sig))
json.dump(res,open("checks.json","w"),indent=1)

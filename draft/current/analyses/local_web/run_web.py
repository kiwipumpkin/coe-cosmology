"""Ensemble driver. usage: python3 run_web.py <s_gpc> <nreal> <seed> <mode: full|tube> [tag]
full: whole region 14.3-30 Gpc (CoE interior excluded).  tube: only cells within T=6s of the CMB ring (holds ~all of the effect, see README).
Saves raw_<tag>.npz (checkpointed every 10 realisations; restarts resume)."""
import sys, os, time, json, numpy as np
from web_core import *
s=float(sys.argv[1]); nreal=int(sys.argv[2]); seed=int(sys.argv[3]); mode=sys.argv[4]
tag=sys.argv[5] if len(sys.argv)>5 else f"s{s:g}_{mode}_seed{seed}"
T=6*s; fn=f"raw_{tag}.npz"
GAP=float(sys.argv[6]) if len(sys.argv)>6 else 0.0   # optional: keep every cell at least GAP Gpc outside the 14.3 Gpc sphere (sensitivity run)
T=6*s+GAP
keys=["f_rmsw","f_left0","f_best0","f_mult0","f_left4","f_best4","f_mult4",
      "u_rmsw","u_left0","u_best0","u_mult0","u_left4","u_best4","u_mult4",
      "n_cells","n_clusters","n_fil","n_wall","med_void_def","med_cluster","med_Rv","min_dist","net_rel"]
res={k:[] for k in keys}; A={"f_amps0":[],"u_amps0":[],"f_amps4":[],"u_amps4":[]}
start=0
if os.path.exists(fn):
    z=np.load(fn); res={k:list(z[k]) for k in keys}; A={k:list(z[k]) for k in A}; start=len(res["f_rmsw"])
t0=time.time()
for i in range(start,nreal):
    rng=np.random.default_rng([seed,i])
    if mode=="full": X,M,rad,info=foam(s,rng); Xu=unstructured(M,rad,rng)
    else: X,M,rad,info=foam_tube(s,rng,T,gap=GAP); Xu=unstructured_tube(M,rad,rng,T,gap=GAP)
    for p,Xp in (("f",X),("u",Xu)):
        w=pattern(Xp,M); res[p+"_rmsw"].append(rms(w)/FLOOR)
        for j,tg in (("0",R),("4",R4)):
            l,b,m=score(w,tg); res[f"{p}_left{j}"].append(l); res[f"{p}_best{j}"].append(b); res[f"{p}_mult{j}"].append(m)
            A[f"{p}_amps{j}"].append(amps(tg-w)/FLOOR)
    res["n_cells"].append(info["n_cells"]); res["n_clusters"].append(info["n_clusters"]); res["n_fil"].append(info["n_fil"]); res["n_wall"].append(info["n_wall"])
    res["med_void_def"].append(info["median_void_deficit_kg"]); res["med_cluster"].append(info["median_cluster_kg"]); res["med_Rv"].append(info["median_void_R_eff_gpc"])
    res["min_dist"].append(info["min_dist_gpc"]); res["net_rel"].append(info["net_mass_rel"])
    if (i+1)%10==0 or i==nreal-1:
        np.savez(fn+".tmp.npz",**{k:np.array(v) for k,v in res.items()},**{k:np.array(v) for k,v in A.items()}); os.replace(fn+".tmp.npz",fn)
        f=np.array(res["f_left0"]); u=np.array(res["u_left0"])
        print(f"{tag} {i+1}/{nreal} t={time.time()-t0:.0f}s foam med={np.median(f):.2f} helped={np.mean(f<BASE0):.3f} rmsw={np.median(res['f_rmsw']):.1f} | unstr med={np.median(u):.2f} helped={np.mean(u<BASE0):.3f} rmsw={np.median(res['u_rmsw']):.1f}",flush=True)

"""Jittered lattice: the best lattice (type, spacing, orientation, phase, amplitude from the search) with every site displaced by an
isotropic Gaussian of RMS 3D length f*a (f = 0, 10, 25, 50%), 1,000 realizations per level. Region, density-weighted masses,
net-zero void rescale and buffer rule are re-applied after the shift (sites can move in/out of the region).
Counts: helped/hurt vs the starting leftover (237.3x alone; 88.8x with the four clumps), and vs the perfect lattice.
run: python3 jitter.py search_alone_with4_buf.json -> jitter.json"""
import json, sys, numpy as np
exec(open("run_search.py").read().split("\nif __name__")[0])
LEVELS_J=[0.0,0.10,0.25,0.50]; NREAL=1000

def jitter_run(best,target):
    kind=best["kind"]; a=best["spacing_gpc"]/RV; e,ph=unpack(kind,best["x"]); A=best["amp"]
    P0,s0=lattice_sites(kind,a,e,ph,clip=False)
    near=np.linalg.norm(P0,axis=1)<=S_MAX+2*a; P0,s0=P0[near],s0[near]
    perfect=best["best"]; base=rms(target)/FLOOR; out={}
    for f in LEVELS_J:
        rng=np.random.default_rng(int(1000*f)+11); placed=[];tuned=[];nsites=[]
        for i in range(NREAL if f>0 else 1):
            P=P0+rng.normal(size=P0.shape)*(f*a/np.sqrt(3)) if f>0 else P0.copy()
            ok=in_region(P); P,s=P[ok],s0[ok]
            for it in range(3):                       # buffer rule at the fitted amplitude (same as the search)
                w,vs,Vs,M=pattern(P,s,a,kind); Ause=float(np.clip(A,-amp_max(vs),amp_max(vs)))
                keep=np.linalg.norm(P,axis=1)>1.0+radius_gpc(np.abs(M*Ause),1.0)/RV
                if keep.all(): break
                P,s=P[keep],s[keep]
            placed.append(rms(target-Ause*w)/FLOOR); tuned.append(best_amp(target,w,vs)[1]/FLOOR); nsites.append(len(P))
        placed=np.array(placed); tuned=np.array(tuned)
        out[f"{int(f*100)}%"]=dict(jitter_rms_gpc=f*best["spacing_gpc"],n=len(placed),
            helped_vs_start=int((placed<base).sum()),hurt_vs_start=int((placed>=base).sum()),
            better_than_perfect=int((placed<perfect-1e-9).sum()),
            median=float(np.median(placed)),p5=float(np.percentile(placed,5)),p95=float(np.percentile(placed,95)),min=float(placed.min()),max=float(placed.max()),
            tuned_helped_vs_start=int((tuned<base).sum()),tuned_median=float(np.median(tuned)),tuned_min=float(tuned.min()),median_sites=float(np.median(nsites)))
        o=out[f"{int(f*100)}%"]; print(kind,best["target"],f"jitter {f:.2f}: helped/hurt vs {base:.1f} = {o['helped_vs_start']}/{o['hurt_vs_start']}  med={o['median']:.1f} min={o['min']:.1f} p95={o['p95']:.1f}  tuned med={o['tuned_median']:.1f}",flush=True)
    return dict(kind=kind,spacing_gpc=best["spacing_gpc"],target=best["target"],perfect=perfect,start=base,levels=out)

if __name__=="__main__":
    res=json.load(open(sys.argv[1])); out={}
    for tname,target in (("alone",R),("with4",R4)):
        rows=[r for r in res if r["target"]==tname]; b=min(rows,key=lambda r:r["best"])
        out[tname]=jitter_run(b,target)
    json.dump(out,open("jitter.json","w"),indent=1)

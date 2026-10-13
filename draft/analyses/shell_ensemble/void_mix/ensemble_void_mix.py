"""with_clumps setup (surround + exp-lock gradient, compensated clumps) with the void share f_v of the random objects varied.
Zero net contrast kept: void deficits rescaled to equal the cluster mass total (same convention as the original 50/50 code).
delta cap: voids start at delta=-0.5 (original); when rescaled they deepen, capped at delta=-1 (empty); beyond that they grow in radius."""
import json, os, sys, time, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); WC=os.path.join(os.path.dirname(HERE),"with_clumps")
g={"__file__":os.path.join(WC,"ensemble_with_clumps.py"),"__name__":"imp"}
src=open(g["__file__"]).read().split('\nif __name__')[0]; os.chdir(WC); exec(src,g); os.chdir(HERE)
field3,amps,rms,RV,FLOOR,KG,LEVELS,radius_gpc,_draw,C,P4,Rcl,baseline,cnt=[g[k] for k in
   ("field3","amps","rms","RV","FLOOR","KG","LEVELS","radius_gpc","_draw","C","P4","Rcl","baseline","cnt")]

def realization(N,lev,rng,sub,fv):
    lo,hi=LEVELS[lev]
    if fv=="rand": fv=rng.uniform(0,1)
    nv=int(np.clip(round(fv*N),1,N-1)); nc=N-nv                        # at least one of each so net contrast can be zero
    Mc=10**rng.uniform(lo,hi,nc); Mv0=10**rng.uniform(lo,hi,nv)
    s=Mc.sum()/Mv0.sum(); Mv=Mv0*s                                       # zero net contrast (same rule as original)
    dv=np.minimum(0.5*s,1.0)                                             # void depth |delta|: 0.5*s, capped at 1 (empty)
    rv=radius_gpc(Mv,dv)                                                 # radius holding the deficit at that depth
    M=np.r_[Mc,-Mv]; rad=np.r_[radius_gpc(Mc,1.0),rv]/RV; rmin=np.ones(N)
    P=np.empty((N,3)); P[:nc]=_draw(nc,rng,C,rmin[:nc],1,True); P[nc:]=_draw(nv,rng,C,rmin[nc:],-1,True)
    if sub:     # redraw random objects that would overlap a placed clump (same rule as with_clumps)
        for _ in range(100):
            dd=np.linalg.norm(P[:,None,:]-P4[sub][None],axis=2); bad=np.where((dd<(rad[:,None]+Rcl[sub][None])).any(1))[0]
            if not bad.size: break
            for kind in (1,-1):
                ii=bad[np.sign(M[bad])==kind]
                if ii.size: P[ii]=_draw(ii.size,rng,C,rmin[ii],kind,True)
    Dd=np.linalg.norm(P,axis=1); w=(M/KG)@field3(np.arctan2(P[:,1],P[:,0]),np.arcsin(P[:,2]/Dd),Dd)
    return w,dict(fv=nv/N,depth=float(dv),max_void_r=float(rv.max()),void_scale=float(s))

def run(N,lev,nreal,Rk,sub,fv,seed):
    rng=np.random.default_rng(seed); res=dict(rms=[],rms_best=[],amps=[],amps_best=[],fv=[],depth=[],max_void_r=[],void_scale=[])
    for i in range(nreal):
        w,info=realization(N,lev,rng,sub,fv); s=float(np.dot(Rk,w)/np.dot(w,w))
        res["rms"].append(rms(Rk-w)); res["rms_best"].append(rms(Rk-s*w)); res["amps"].append(amps(Rk-w)); res["amps_best"].append(amps(Rk-s*w))
        for k,v in info.items(): res[k].append(v)
    return {k:np.array(v) for k,v in res.items()}

FVS=[0.10,0.25,0.50,0.75,0.90,"rand"]
if __name__=="__main__":
    NREAL=int(sys.argv[1]); ks=[int(x) for x in sys.argv[2].split(",")]; Ns=[int(x) for x in sys.argv[3].split(",")]; levs=sys.argv[4]
    tag=sys.argv[5]; out={}; t0=time.time()
    for k in ks:
        Rk,sub=baseline(k); base=rms(Rk)/FLOOR; A0=amps(Rk)
        for N in Ns:
            for lev in levs:
                for j,fv in enumerate(FVS):
                    seed=1+3*(N==500)+"abc".index(lev)   # same seeds as with_clumps (fv=50% reproduces it; common random numbers across f_v)
                    res=run(N,lev,NREAL,Rk,sub,fv,seed); fl="rand" if fv=="rand" else f"{int(round(fv*100))}"
                    key=f"k{k}_N{N}_{lev}_fv{fl}"; np.savez_compressed(f"raw_{key}.npz",**res)
                    out[key]=dict(k=k,N=N,level=lev,fv=fl,subset=sub,baseline=base,placed=cnt(res["rms"],res["amps"],base,A0),best=cnt(res["rms_best"],res["amps_best"],base,A0),
                                  void_depth_median=float(np.median(res["depth"])),void_depth_capped_frac=float(np.mean(res["depth"]>=1.0)),
                                  max_void_radius_gpc=float(res["max_void_r"].max()),median_void_scale=float(np.median(res["void_scale"])))
                    p=out[key]["placed"]; b=out[key]["best"]
                    print(key,f"base={base:.1f} help/hurt={p['rms_lower']}/{p['rms_higher']} med={p['median']:.1f} min={p['min']:.1f} best_med={b['median']:.1f} best_min={b['min']:.1f} "
                          f"quad={p['m35_lower_and_m12_not_higher']} strict={p['each_m3_m4_m5_lower_and_each_m1_m2_not_higher']} depth={out[key]['void_depth_median']:.2f} capped={out[key]['void_depth_capped_frac']:.2f} rmax={out[key]['max_void_radius_gpc']:.2f} t={time.time()-t0:.0f}s",flush=True)
    json.dump(out,open(f"part_{tag}.json","w"),indent=1)

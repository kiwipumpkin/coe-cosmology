"""Surround+gradient random web (identical to ../surround grad_bare) on top of baselines with k=1..4 fitted clumps (compensated)."""
import json, os, sys, time, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__))
g={"__file__":os.path.join(os.path.dirname(HERE),"surround","ensemble_surround.py"),"__name__":"imp"}
exec(open(g["__file__"]).read().split('\nif __name__')[0],g); os.chdir(HERE)
field3,amps,R,rms,RV,FLOOR,KG,LEVELS,radius_gpc,_draw,coe_vec,phi,B_SH=[g[k] for k in
   ("field3","amps","R","rms","RV","FLOOR","KG","LEVELS","radius_gpc","_draw","coe_vec","phi","B_SH")]
B=np.load("baseline_patterns.npz"); F4,kk,Fc,P4=B["F4"],B["kk"],B["Fcomp_unit"],B["P4"]
bj=json.load(open("baselines.json")); ELL=phi+180.0; C=coe_vec(ELL)
Rcl=radius_gpc(kk*KG,1.0)/RV        # delta=+1 radius of each placed clump (~2.8 Gpc)

def baseline(k):
    if k==0: return R.copy(),[]
    sub=bj["baselines"][str(k)]["best"]["sub"]
    return R-(kk[sub]@F4[sub]-kk[sub].sum()*Fc), sub

def realization(N,lev,rng,sub):
    lo,hi=LEVELS[lev]; nh=N//2
    Mc=10**rng.uniform(lo,hi,nh); Mv=10**rng.uniform(lo,hi,N-nh); Mv*=Mc.sum()/Mv.sum(); M=np.r_[Mc,-Mv]
    rad=np.r_[radius_gpc(Mc,1.0),radius_gpc(Mv,0.5)]/RV; rmin=np.ones(N)
    P=np.empty((N,3)); P[:nh]=_draw(nh,rng,C,rmin[:nh],1,True); P[nh:]=_draw(N-nh,rng,C,rmin[nh:],-1,True)
    nred=0
    if sub:     # exclude overlaps with placed clumps: redraw any object whose sphere touches a clump's sphere
        for _ in range(100):
            dd=np.linalg.norm(P[:,None,:]-P4[sub][None],axis=2); bad=np.where((dd<(rad[:,None]+Rcl[sub][None])).any(1))[0]
            if not bad.size: break
            nred+=bad.size
            for kind in (1,-1):
                ii=bad[np.sign(M[bad])==kind]
                if ii.size: P[ii]=_draw(ii.size,rng,C,rmin[ii],kind,True)
    Dd=np.linalg.norm(P,axis=1); w=(M/KG)@field3(np.arctan2(P[:,1],P[:,0]),np.arcsin(P[:,2]/Dd),Dd)
    return w,M,P,nred

def run(N,lev,nreal,Rk,sub,seed):
    rng=np.random.default_rng(seed); res=dict(rms=[],rms_best=[],amps=[],amps_best=[],nred=[]); keep=None
    for i in range(nreal):
        w,M,P,nred=realization(N,lev,rng,sub)
        s=float(np.dot(Rk,w)/np.dot(w,w)); res["rms"].append(rms(Rk-w)); res["rms_best"].append(rms(Rk-s*w))
        res["amps"].append(amps(Rk-w)); res["amps_best"].append(amps(Rk-s*w)); res["nred"].append(nred)
        if i==0: keep=dict(M=M.tolist(),P=(P*RV).tolist())
    return {k:np.array(v) for k,v in res.items()},keep

def cnt(r,A,base,A0):
    r=r/FLOOR; q35=np.linalg.norm(A[:,2:5],axis=1); q12=np.linalg.norm(A[:,0:2],axis=1); q35_0=np.linalg.norm(A0[2:5]); q12_0=np.linalg.norm(A0[0:2])
    lo=int((r<base).sum())
    return dict(n=len(r),rms_lower=lo,rms_higher=len(r)-lo,m35_lower=int((q35<q35_0).sum()),m12_higher=int((q12>q12_0).sum()),
        m35_lower_and_m12_not_higher=int(((q35<q35_0)&(q12<=q12_0)).sum()),
        each_m3_m4_m5_lower_and_each_m1_m2_not_higher=int((np.all(A[:,2:5]<A0[2:5],1)&np.all(A[:,0:2]<=A0[0:2],1)).sum()),
        median=float(np.median(r)),p5=float(np.percentile(r,5)),p95=float(np.percentile(r,95)),min=float(r.min()),max=float(r.max()))

if __name__=="__main__":
    NREAL=int(sys.argv[1]) if len(sys.argv)>1 else 200
    out=dict(baselines={},counts={},examples={}); t0=time.time()
    for k in (0,1,2,3,4):
        Rk,sub=baseline(k); base=rms(Rk)/FLOOR; A0=amps(Rk)
        out["baselines"][k]=dict(subset=sub,rms_over_floor=base,amps_over_floor=(A0/FLOOR).tolist())
        seed=1
        for N in (200,500):
            for lev in "abc":
                res,keep=run(N,lev,NREAL,Rk,sub,seed); seed+=1     # same seeds as ../surround grad_bare (k=0 reproduces it)
                key=f"k{k}_N{N}_{lev}"; np.savez_compressed(f"raw_{key}.npz",**res)
                out["counts"][key]=dict(placed=cnt(res["rms"],res["amps"],base,A0),best=cnt(res["rms_best"],res["amps_best"],base,A0),
                                        mean_redrawn=float(res["nred"].mean()))
                if lev=="b" and N==200: out["examples"][key]=keep
                p=out["counts"][key]["placed"]; b=out["counts"][key]["best"]
                print(key,f"base={base:.1f} lower/higher={p['rms_lower']}/{p['rms_higher']} med={p['median']:.1f} min={p['min']:.1f} best_med={b['median']:.1f} best_min={b['min']:.1f} quad={p['m35_lower_and_m12_not_higher']} strict={p['each_m3_m4_m5_lower_and_each_m1_m2_not_higher']} redraw={res['nred'].mean():.2f} t={time.time()-t0:.0f}s",flush=True)
    json.dump(out,open("counts.json","w"),indent=1)

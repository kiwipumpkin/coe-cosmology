"""Random web (void_mix setup: surround + exp-lock gradient, zero net contrast, overlap redraw) on top of the refit baselines:
masses only vs masses+voids, k=1..4, N=200, level b, void share 50% and 90%, 1,000 layouts, seed as void_mix (common random numbers)."""
import json, os, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); VM=os.path.dirname(HERE)
vm={"__file__":os.path.join(VM,"ensemble_void_mix.py"),"__name__":"imp"}
src=open(vm["__file__"]).read().split('\nif __name__')[0]; os.chdir(VM); exec(src,vm); os.chdir(HERE)
fm={"__file__":os.path.join(HERE,"fit_masses_voids.py"),"__name__":"imp"}
exec(open(fm["__file__"]).read().split('\nif __name__')[0],fm); os.chdir(HERE)
rms,amps,FLOOR,RV,cnt,R=[vm[k] for k in ("rms","amps","FLOOR","RV","cnt","g")][:5]+[fm["R"]]
C=vm["C"]
F=json.load(open("fit_masses_voids.json")); out={}
for k in (1,2,3,4):
    for var in ("masses_only","masses_voids"):
        r=F[f"k{k}_{var}"]; q=np.array(r["q"]); signed=(var=="masses_voids")
        Rk=R-fm["model"](q,k,signed); base=rms(Rk)/FLOOR; A0=amps(Rk)
        a=np.radians(r["ell_deg"]); D=np.array(r["dist_gpc"])/RV
        vm["P4"]=np.c_[D*np.cos(a),D*np.sin(a),np.zeros(k)]; vm["Rcl"]=np.array(r["radius_gpc_at_abs_delta_1"])/RV
        rC=np.linalg.norm(vm["P4"]-C,axis=1)*RV
        for fv in (0.50,0.90):
            res=vm["run"](200,"b",1000,Rk,list(range(k)),fv,2)
            key=f"k{k}_{var}_fv{int(fv*100)}"; np.savez_compressed(f"raw_{key}.npz",**res)
            out[key]=dict(k=k,variant=var,fv=int(fv*100),baseline=base,dist_from_CoE_gpc=rC.tolist(),
                          placed=cnt(res["rms"],res["amps"],base,A0),best=cnt(res["rms_best"],res["amps_best"],base,A0))
            p=out[key]["placed"]; b=out[key]["best"]
            print(key,f"base={base:.1f} help/hurt={p['rms_lower']}/{p['rms_higher']} med={p['median']:.1f} min={p['min']:.1f} best_med={b['median']:.1f} best_min={b['min']:.1f} quad={p['m35_lower_and_m12_not_higher']} strict={p['each_m3_m4_m5_lower_and_each_m1_m2_not_higher']} rCoE={np.round(rC,1).tolist()}",flush=True)
json.dump(out,open("ensemble_fitted.json","w"),indent=1)

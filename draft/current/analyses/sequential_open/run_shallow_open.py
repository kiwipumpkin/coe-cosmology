"""Shallow-void runs with the open region: 'grow cmb_aware <seed>' (B) or 'near <nseeds>' (A). Original shallow code, region patched.
B cap: original 6000 voids scaled by the region volume (x2.364 -> 14,184). Checkpoints."""
import os,sys,time,json,numpy as np
from common import *
SH=f"{SEQ}/shallow"; sys.path.insert(0,SH); os.chdir(SH)
import shallow_core as sc
os.chdir(OUT)
patch(sc)
from shallow_core import bank_patterns, def_mass, Packer, rand_dir, check_no_overlap, ldiff, vol, LOCAL_RMAX, S_MAX_G
CAP=int(round(6000*V_REG/V_OLD_B))
mode=sys.argv[1]; os.makedirs(f"{OUT}/ckpt",exist_ok=True)
if mode=="grow":
    rule,seed=sys.argv[2],int(sys.argv[3]); Rv=1.5; CK=f"{OUT}/ckpt/shallow_{rule}_s{seed}{SUF}.pkl"
    if os.path.exists(CK): st=pickle.load(open(CK,"rb")); print("resumed",len(st["order"]),flush=True)
    else: st=dict(rng=np.random.default_rng(seed),pk=Packer(h=4.0),fronts=[],order=[],w=np.zeros(360),fail=0)
    rng,pk,fronts,order,w,fail=st["rng"],st["pk"],st["fronts"],st["order"],st["w"],st["fail"]; t0=time.time(); last=time.time()
    while fail<3 and len(order)<CAP:
        cands=[]
        for _ in range(600):
            if rule!="random" and fronts:
                j=rng.integers(len(fronts)); x=fronts[j]+rand_dir(rng)*2*Rv*rng.uniform(1.02,1.12)
            else: x=O.rand_in_region(rng,RMIN,S_MAX_G,1)[0]
            if O.region_ok(x,Rv)[0] and pk.free(x,Rv): cands.append(x)
            if len(cands)>=(30 if rule=="cmb_aware" else 1): break
        if not cands: fronts=fronts[len(fronts)//2:]; fail+=1; continue
        fail=0
        if rule=="cmb_aware":
            F=bank_patterns(np.array(cands),Rv); i=int(np.argmin((((R-w)[None]-F)**2).sum(1))); x=cands[i]; w=w+F[i]
        else: x=cands[0]
        pk.add(x,Rv,0.0,"void"); fronts.append(x); order.append(x)
        if len(order)%200==0 and time.time()-last>60:
            save_atomic(dict(rng=rng,pk=pk,fronts=fronts,order=order,w=w,fail=fail),CK); last=time.time()
            print(rule,seed,len(order),f"{time.time()-t0:.0f}s ckpt",flush=True)
    save_atomic(dict(rng=rng,pk=pk,fronts=fronts,order=order,w=w,fail=fail,done=True),CK)
    X=np.array(order); F=bank_patterns(X,Rv); C=np.cumsum(F,0); left=np.sqrt(np.mean((R[None]-C)**2,1))/FLOOR; wf=C[-1]
    Xa,ra,_,_=pk.arrays()
    out=dict(rule=rule,seed=seed,n_void=len(X),cap=CAP,hit_cap=len(X)>=CAP,void_vol_frac=float(len(X)*vol(Rv)/V_REG),final=float(left[-1]),min=float(left.min()),
             argmin=int(np.argmin(left))+1,min_gap=float(check_no_overlap(Xa,ra,0.0)),where_voids=where(X),final_amps=(amps(R-wf)/FLOOR)[:6].tolist(),
             traj_obj=list(range(1,len(left)+1,5)),traj_left=left[::5].round(3).tolist())
    json.dump(out,open(f"{OUT}/shallow_{rule}_s{seed}{SUF}.json","w")); np.savez(f"{OUT}/final_shallow_{rule}_s{seed}{SUF}.npz",X=X,Rv=Rv,w=wf,left=left[-1])
    print("DONE",{q:v for q,v in out.items() if not q.startswith("traj")},flush=True)
else:
    # A near fix, shallow voids: original near() with banks scaled to keep candidate density (x5.16)
    src=open(f"{SH}/run_shallow.py").read().split('if __name__')[0]
    g={"__file__":f"{SH}/run_shallow.py"}; os.chdir(SH); exec(src,g); os.chdir(OUT)
    g["region_ok"]=O.region_ok; g["rand_in_shell"]=O.rand_in_region
    nb_=int(round(5000*V_NEW_A/V_OLD_A)); ns=int(sys.argv[2]); runs=[]
    for s in range(ns):
        fn=f"{OUT}/ckpt/shallowA_seed{s}.json"
        if os.path.exists(fn): r=json.load(open(fn))
        else: r=g["near"](s,nbank=nb_); json.dump(r,open(fn,"w"))
        runs.append(r); print("seed",s,"left",round(r["left"],2),"voids",r["n"],"where",where(np.array(r["X"])),flush=True)
    b=min(runs,key=lambda r:r["left"])
    json.dump(dict(best=b,nbank=nb_,runs=[{k:v for k,v in r.items() if k!="w"} for r in runs]),open(f"{OUT}/shallowA_near.json","w"))
    print("DONE A shallow best",b["seed"],b["left"],flush=True)

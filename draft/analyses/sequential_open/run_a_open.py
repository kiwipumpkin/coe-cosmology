"""Scenario A near fix (deep voids + clusters) with the open region. Original local_fit, region patched, candidate bank x5.16."""
import os,sys,json,numpy as np
from common import *
os.chdir(SEQ); import scenario_a as sa; os.chdir(OUT)
patch(sa)
nb_=int(round(6000*V_NEW_A/V_OLD_A)); ns=int(sys.argv[1]) if len(sys.argv)>1 else 4; runs=[]
os.makedirs(f"{OUT}/ckpt",exist_ok=True)
for sd in range(ns):
    fn=f"{OUT}/ckpt/deepA_seed{sd}.npz"
    if os.path.exists(fn):
        z=np.load(fn,allow_pickle=True); r=dict(seed=sd,X=z["X"],r=z["r"],M=z["M"],kind=z["kind"],w=z["w"],left=float(z["left"]),traj=json.loads(str(z["traj"])),margin=float(z["margin"]))
    else:
        r=sa.local_fit(seed=sd,nbank=nb_)
        np.savez(fn,X=r["X"],r=r["r"],M=r["M"],kind=r["kind"],w=r["w"],left=r["left"],traj=json.dumps(r["traj"]),margin=r["margin"])
    runs.append(r); print(f"seed {sd}: {r['left']:.2f}x objs={len(r['M'])} gap={r['margin']:.3f} where={where(r['X'])}",flush=True)
b=min(runs,key=lambda r:r["left"])
np.savez(f"{OUT}/deepA_best.npz",X=b["X"],r=b["r"],M=b["M"],kind=b["kind"],w=b["w"],left=b["left"])
json.dump(dict(best_seed=b["seed"],nbank=nb_,runs=[dict(seed=r["seed"],left=r["left"],n_obj=len(r["M"]),min_gap=r["margin"],traj=[t["left"] for t in r["traj"]],
     where=where(r["X"]),final_amps=(amps(R-r["w"])/FLOOR)[:6].tolist()) for r in runs]),open(f"{OUT}/deepA_near.json","w"),indent=1)
print("DONE A deep best",b["seed"],b["left"],flush=True)

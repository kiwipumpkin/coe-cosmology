"""Scenario B (deep voids + clusters) with the open region. Original Grow code, region patched; checkpoints every 100 voids.
usage: python3 run_b_open.py <rule> <seed>   -> ckpt/b_<rule>_s<seed>.pkl (resume), b_<rule>_s<seed>.json, final_b_<rule>_s<seed>.npz"""
import os,sys,time,json,numpy as np
from common import *
os.chdir(SEQ); import scenario_b as sb; os.chdir(OUT)
patch(sb)
rule,seed=sys.argv[1],int(sys.argv[2]); Rv=1.5
os.makedirs(f"{OUT}/ckpt",exist_ok=True); CK=f"{OUT}/ckpt/b_{rule}_s{seed}{SUF}.pkl"
if os.path.exists(CK):
    st=pickle.load(open(CK,"rb")); g=st["g"]; fail=st["fail"]; print("resumed",len(g.voids),"voids",len(g.order),"objs",flush=True)
else:
    g=sb.Grow(rule,seed,Rv); fail=0
t0=time.time(); last=time.time()
s=g
while fail<20:
    r=s.step_void()
    if r is None:
        if s.fronts: s.fronts=s.fronts[len(s.fronts)//2:]
        fail+=1; continue
    x,Mv=r; fail=0
    Rc=sb.clu_R(-Mv/sb.NRIM,x); n=0
    for _ in range(sb.NRIM):
        if s.step_cluster(x,-Mv/sb.NRIM,Rc): n+=1
    if n<sb.NRIM: s.undo_unit(1+n); continue
    if len(s.voids)%100==0 and time.time()-last>60:
        save_atomic(dict(g=s,fail=fail),CK); last=time.time()
        print(f"{rule} s{seed} voids={len(s.voids)} objs={len(s.order)} t={time.time()-t0:.0f}s ckpt",flush=True)
save_atomic(dict(g=s,fail=fail,done=True),CK)
X,M,k,left,w=sb.trajectory(s.order); Xa,ra,Ma,ka=s.pk.arrays()
gap=seq_core.check_no_overlap(Xa,ra,0.0)
idx=np.r_[np.arange(min(2000,len(M))),np.arange(2000,len(M),10)].astype(int) if len(M)>2000 else np.arange(len(M))
nv=np.cumsum(k=="void"); vf=float(np.sum(k=="void")*seq_core.vol(Rv)/V_REG)
iv=np.array(ka)=="void"
out=dict(rule=rule,seed=seed,Rv=Rv,n_obj=int(len(M)),n_void=int((k=="void").sum()),n_clu=int((k=="clu").sum()),void_vol_frac=vf,net_kg=float(M.sum()),
         min_gap_gpc=float(gap),final_left=float(left[-1]),min_left=float(left.min()),argmin=int(np.argmin(left))+1,web_rms=float(rms(w)/FLOOR),
         final_amps=(amps(R-w)/FLOOR)[:6].tolist(),where_all=where(Xa),where_voids=where(Xa[iv]),where_clusters=where(Xa[~iv]),
         traj_obj=(idx+1).tolist(),traj_left=left[idx].round(3).tolist())
json.dump(out,open(f"{OUT}/b_{rule}_s{seed}{SUF}.json","w"))
np.savez(f"{OUT}/final_b_{rule}_s{seed}{SUF}.npz",Xp=Xa,rp=ra,kp=np.array(ka).astype(str),Mp=Ma,w=w,left=left[-1],X=X,M=M,kind=k.astype(str))
print("DONE",{q:v for q,v in out.items() if not q.startswith("traj")},flush=True)

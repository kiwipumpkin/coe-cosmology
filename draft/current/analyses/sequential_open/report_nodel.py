"""Report for the no-deletion minimal edit (minedit2.py NODEL=1): reloc = near-ring bodies moved to the nearest legal spots first,
forced = near-ring bodies fitted out of the 3 Gpc band from their own spots; then minimal move/resize of other bodies."""
import os,json,pickle,glob,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import *
from check_layout import check
from edit_before_after import figure as before_after, run_figure
plt.rcParams.update({"font.size":7.5})
WEBS=[("b_random_s1_tidalcore","random s1","#95a5a6"),("b_random_s2_tidalcore","random s2","#717d7e"),("b_random_s3_tidalcore","random s3","#2c3e50"),
      ("b_foam_rim_s1_tidalcore","foam s1","#85c1e9"),("b_foam_rim_s2_tidalcore","foam s2","#3498db"),("b_foam_rim_s3_tidalcore","foam s3","#1b4f72")]
MODES=[("reloc","relocate first","-"),("forced","forced out (fitted)","--")]
D={}
for t,nm,c in WEBS:
    for m,ml,ls in MODES:
        k=f"{t}_nodel_{m}"; fj=f"{OUT}/minedit2/{k}.json"; ck=f"{OUT}/ckpt/minedit2_{k}.pkl"
        if os.path.exists(fj): D[k]=json.load(open(fj)); D[k]["final"]=True
        elif os.path.exists(ck): D[k]=dict(curve=pickle.load(open(ck,"rb"))["curve"],final=False)
fig,ax=plt.subplots(figsize=(7,6.8))
for t,nm,c in WEBS:
    for m,ml,ls in MODES:
        k=f"{t}_nodel_{m}"
        if k in D: cv=np.array(D[k]["curve"]); ax.plot(np.maximum(cv[:,0],0.8),cv[:,1],ls,color=c,lw=1.4,label=f"{nm}, {ml}"+("" if D[k]["final"] else " (running)"))
ax.axhline(1,color="green",lw=1.5); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("existing bodies moved / resized (log; first jump = all near-ring bodies moved out)"); ax.set_ylabel("leftover (× 1×)")
ax.set_title("Minimal edit, no deletion: every body kept, near-ring bodies moved ≥3 Gpc from the ring",fontsize=8.5)
hd=[l for l in ax.get_lines() if not l.get_label().startswith("_")]+[Line2D([],[],color="green",lw=1.5,label="1× (CMB level)")]
fig.legend(handles=hd,loc="upper center",ncol=2,frameon=False,fontsize=6.5,bbox_to_anchor=(0.5,0.0)); fig.tight_layout(); fig.savefig(f"{OUT}/nodel_curve.png",dpi=200,bbox_inches="tight"); plt.close(fig)
S={}
for k,d in D.items():
    if not d["final"]: S[k]=dict(running=True,last=d["curve"][-1]); continue
    mo=np.array(d["moved"],bool); rs=np.array(d["resized"],bool); nr=np.array(d["near_ring"],bool); mv=np.array(d["move"]); a=d["after"]; ed=mo|rs
    S[k]=dict(n_bodies=d["n_bodies"],edited=int(ed.sum()),frac=round(float(ed.sum()/d["n_bodies"]),4),moved=int(mo.sum()),resized=int(rs.sum()),near_ring=d["n_near_ring"],
      before=round(d["before"]["total"],1),left=round(d["left"],3),reached=d["reached"],m1_12=round(float(np.hypot(a["rms_m1_5"],a["rms_m6_12"])),3),m13=round(a["rms_13up"],3),sym=round(a["sym"],2),
      move_near_ring_median=round(float(np.median(mv[nr&mo])),2) if (nr&mo).any() else None,move_near_ring_p90=round(float(np.percentile(mv[nr&mo],90)),2) if (nr&mo).any() else None,
      move_other_median=round(float(np.median(mv[~nr&mo])),2) if (~nr&mo).any() else None,move_other_p90=round(float(np.percentile(mv[~nr&mo],90)),2) if (~nr&mo).any() else None,
      mass_scale_10_90=[round(float(x),2) for x in np.percentile(np.array(d["a_new"])[ed],[10,90])],
      dist_us_median=round(float(np.median(np.array(d["dist_us_new"])[ed])),1),dist_ring_median=round(float(np.median(np.array(d["dist_ring_new"])[ed])),1),
      net_mass_change=d["net_mass_change"],legality=check(f"{OUT}/minedit2/{k}.json",k.split("_nodel")[0]))
json.dump(S,open(f"{OUT}/minedit2/summary_nodel.json","w"),indent=1)
for k,v in S.items(): print(k,v)
fin=[k for k in D if D[k]["final"]]
if fin:
    pick=[]
    for fam in ("random","foam"):
        for m in ("reloc","forced"):
            c=[k for k in fin if fam in k and k.endswith(m)]
            if c: pick.append(min(c,key=lambda k:(not D[k]["reached"],D[k]["n_edited"] if D[k]["reached"] else D[k]["left"])))
    col={t:c for t,_,c in WEBS}; nmd={t:nm for t,nm,_ in WEBS}
    items=[(f"{OUT}/minedit2/{k}.json",f"{nmd[k.split('_nodel')[0]]}, no deletion, {dict((m,ml) for m,ml,_ in MODES)[k.split('_nodel_')[1]]}",col[k.split('_nodel')[0]]) for k in pick]
    before_after(items,f"{OUT}/nodel_side_views.png")

# default (10 Oct): one image per run
for t,nm,c in WEBS:
    for m,ml,_ in MODES:
        k=f"{t}_nodel_{m}"
        if k in D and D[k]["final"]:
            short=t.replace("b_","").replace("_rim","").replace("_tidalcore","")
            run_figure(f"{OUT}/minedit2/{k}.json",f"{nm}, no deletion, {ml}",c,f"{OUT}/nodel_run_{short}_{m}.png")

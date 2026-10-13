"""Report for the entropy-constrained minimal edit (minedit_ent.py): table, one PNG per run, ring-profile chart for the best run."""
import json,os,glob,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import edit_before_after as EB
D="minedit_ent"; runs=[]
for t in ("b_random_s1","b_random_s2","b_random_s3","b_foam_rim_s1","b_foam_rim_s2","b_foam_rim_s3"):
    p=f"{D}/{t}_tidalcore_entropy.json"
    if os.path.exists(p): runs.append((t,json.load(open(p)),p))
COL={"random":"#2c3e50","foam":"#8e44ad"}
lab=lambda t: ("random" if "random" in t else "foam")+" growth, seed "+t[-1]
K=("radial","ring","angular","link_cc","link_cv","link_vv","ncount")
rows=["| web | edits (% of web) | moved | resized | leftover before → after | ≤1x? | after m1–12 / m13+ | after symmetric | entropy all ≥ start? | within 3 Gpc of ring start→end |","|---|---|---|---|---|---|---|---|---|---|"]
erows=["| web | "+" | ".join(K)+" |","|---|"+"---|"*len(K)]
for t,d,p in runs:
    sa=d["after"]; m112=np.hypot(sa["rms_m1_5"],sa["rms_m6_12"])
    rows.append(f"| {lab(t)} | {d['n_edited']} ({100*d['frac_edited']:.1f}%) | {d['n_moved']} | {d['n_resized']} | {d['before']['total']:.0f}× → {d['left']:.2f}× | {'yes' if d['left']<=1 else 'no'} | {m112:.2f}× / {sa['rms_13up']:.2f}× | {100*sa['sym']:.0f}% | {'yes' if d['entropy_ok'] else 'NO'} | {d['n_start_within3_ring']}→{d['n_end_within3_ring']} |")
    erows.append(f"| {lab(t)} | "+" | ".join(f"{d['entropy_before'][k]:.4f}→{d['entropy_after'][k]:.4f}" for k in K)+" |")
    EB.run_figure(p,lab(t)+" — moves/resizes under the entropy rule (no ring backstop)",COL["random" if "random" in t else "foam"],f"entropy_run_{t}.png")
open(f"{D}/table.md","w").write("\n".join(rows)+"\n\nEntropy before → after (radial/ring/angular: volume-weighted, 0 = uniform in volume; others: Shannon of the histogram)\n\n"+"\n".join(erows)+"\n")
print("\n".join(rows)); print("\n".join(erows))
if runs:
    t,d,_=min(runs,key=lambda r:r[1]["left"]); rp=d["ring_profile"]; e=np.array(rp["edges"]); c=0.5*(e[1:]+e[:-1]); v=np.array(rp["vol"])
    fig,ax=plt.subplots(figsize=(8,5.6))
    ax.plot(c,np.array(rp["before"])/v,"o-",color="0.55",lw=2,label="before (as grown)")
    ax.plot(c,np.array(rp["after"])/v,"s--",color="#1e8449",lw=2,label=f"after ({d['n_edited']} edits)")
    ax.set_xlabel("distance of body centre from the CMB ring (Gpc)",fontsize=12); ax.set_ylabel("bodies per Gpc³",fontsize=12)
    ax.axhline(d["n_bodies"]/(4/3*np.pi*30**3),color="#c0392b",ls=":",lw=1.8,label="uniform-in-volume average")
    ax.set_title(f"Bodies vs distance from the CMB ring: {lab(t)} (best run)",fontsize=13); ax.tick_params(labelsize=11); ax.set_ylim(0,None); ax.grid(alpha=0.3)
    ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.14),ncol=3,frameon=False,fontsize=11)
    fig.savefig("entropy_ring_density.png",dpi=130,bbox_inches="tight"); print("best",t)

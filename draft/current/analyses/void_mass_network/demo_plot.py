import os,sys,glob,json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ME=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(ME,"..","joint_network"))
_c=os.getcwd(); import jn_core as J; os.chdir(_c)
runs={}
for f in glob.glob(os.path.join(ME,"results","demo_*.json")):
    for r in json.load(open(f)): runs.setdefault((r["Mmin"],r["mode"]),[]).append(r)
lab={(3e15,"lcdm"):("clusters >3e15 Msun, cells compensated to LCDM σ (default)","C0"),
     (2e15,"lcdm"):("clusters >2e15 Msun, cells compensated to LCDM σ","C2"),
     (3e15,"exact"):("clusters >3e15 Msun, every cell exactly compensated","C1"),
     (3e15,"hamaus"):("clusters >3e15 Msun, Hamaus+14 trend extrapolated (big voids under-compensated)","C3")}
summ={}
fig,ax=plt.subplots(2,1,figsize=(11,12.5),gridspec_kw=dict(height_ratios=[1.5,1]))
a=ax[0]; bins=np.logspace(-0.3,3.3,90)
for k,(t,c) in lab.items():
    if k not in runs: continue
    L=np.array([r["left"] for r in runs[k]]); P=np.array([r["pattern_rms"] for r in runs[k]])
    summ[f"{k[0]:.0e}_{k[1]}"]=dict(n=len(L),left_median=float(np.median(L)),left_min=float(L.min()),left_max=float(L.max()),
        pattern_rms_median=float(np.median(P)),inside_median=float(np.median([r["pattern_rms_inside"] for r in runs[k]])),
        outside_median=float(np.median([r["pattern_rms_outside"] for r in runs[k]])),frac_below_237=float((L<J.rms(J.R)/J.FLOOR).mean()))
    a.hist(L,bins,color=c,alpha=0.55,label=f"{t}  (n={len(L)}, median {np.median(L):.0f}×)")
a.axvline(1,color="green",lw=2.5,label="1× target (CMB flat)"); a.axvline(237.26,color="k",lw=2,ls="--",label="237× = leftover with no structure")
a.set_xscale("log"); a.set_xlim(0.5,2000); a.set_xlabel("leftover after adding the random network (× CMB floor)"); a.set_ylabel("number of random networks")
a.set_title("Demo null ensemble: rule-following random void–cluster networks through the CMB kernel",fontsize=13)
a.legend(loc="upper center",bbox_to_anchor=(0.5,-0.13),fontsize=10,frameon=False)
# ring panel: residual (target minus network pattern) for the saved default example nearest the median
b=ax[1]; k=(3e15,"lcdm"); deg=np.degrees(J.TH)
if k in runs:
    med=np.median([r["left"] for r in runs[k]]); best=None
    for r in runs[k]:
        fn=os.path.join(ME,"results",f"pattern_{r['Mmin']:.0e}_{r['mode']}_{r['seed']}.npy")
        if os.path.exists(fn) and (best is None or abs(r["left"]-med)<abs(best[0]["left"]-med)): best=(r,np.load(fn))
    if best:
        r,P=best; b.plot(deg,J.R/J.FLOOR,color="0.6",lw=1.2,label="CMB residual to cancel (no structure, 237×)")
        b.plot(deg,(J.R-P)/J.FLOOR,color="#34495e",lw=1.6,label=f"leftover with random network #{r['seed']} ({r['left']:.0f}×)")
b.axhspan(-1,1,color="green",alpha=0.22,lw=0,label="±1× CMB band (thin at this scale)"); b.axhline(0,color="0.5",lw=0.8)
b.set_xlim(0,360); b.set_xticks([0,90,180,270,360]); b.set_xlabel("position around the ring (deg)"); b.set_ylabel("leftover (× floor)")
b.set_title("Around the ring: one typical random network",fontsize=12)
b.legend(loc="upper center",bbox_to_anchor=(0.5,-0.2),fontsize=10,frameon=False,ncol=1)
fig.tight_layout(); fig.savefig(os.path.join(ME,"demo_leftover_hist.png"),dpi=110)
json.dump(summ,open(os.path.join(ME,"demo_summary.json"),"w"),indent=1); print(json.dumps(summ,indent=1))

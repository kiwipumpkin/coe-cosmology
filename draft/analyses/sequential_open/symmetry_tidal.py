"""Mirror symmetry of the combined ring residual (R - bodies) about the CoE axis l = 15.16 / 195.16, modes m = 1-12.
symmetric part = cos(m(l - l0)) terms, antisymmetric = sin terms. sym fraction = sum a^2 / sum (a^2 + b^2)."""
import os,json,glob,re,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import *
AX=np.radians(15.16); TH=np.radians(np.arange(360)); m=np.arange(1,13)
def symfrac(v):
    a=2*np.mean(v[None]*np.cos(m[:,None]*(TH-AX)),1); b=2*np.mean(v[None]*np.sin(m[:,None]*(TH-AX)),1)
    return a,b
a0,b0=symfrac(R); print("start residual: sym fraction %.6f"%((a0**2).sum()/((a0**2+b0**2).sum())))
SC={"b_random":("B random","#7f8c8d"),"b_foam_rim":("B foam rule","#5dade2"),"b_cmb_aware":("B CMB-aware","#7d3c98"),"shallow_cmb_aware":("B CMB-aware, shallow","#b48fc9")}
VAR={"tidal":("our sphere + shell to 30 Gpc filled, core empty","o")}
rows=[]
for fn in glob.glob(f"{OUT}/final_*_s*.npz"):
    b=os.path.basename(fn)[6:-4]; mm=re.match(r"(b_random|b_foam_rim|b_cmb_aware|shallow_cmb_aware)_s(\d+)(_tidal)$",b)
    if not mm: continue
    sc,seed,var=mm.group(1),int(mm.group(2)),(mm.group(3) or "_both")[1:]
    f=np.load(fn); w=f["w"]; res=R-w; a,bb=symfrac(res); aw,bw=symfrac(w)
    P=a**2+bb**2; rows.append(dict(scen=sc,var=var,seed=seed,left=float(rms(res)/FLOOR),sym=float((a**2).sum()/P.sum()),
        sym_m=(a**2/np.maximum(P,1e-300)).tolist(),amp_m=(np.sqrt(P)/FLOOR).tolist(),antisym_rms=float(np.sqrt((bb**2).sum()/2)/FLOOR),sym_rms=float(np.sqrt((a**2).sum()/2)/FLOOR),
        body_sym=float((aw**2).sum()/(aw**2+bw**2).sum())))
json.dump(rows,open(f"{OUT}/symmetry_tidal.json","w"),indent=1)
# chart 1: sym fraction vs leftover, every final layout
fig,ax=plt.subplots(figsize=(7.1,5.6))
for r in rows: ax.scatter(r["left"],r["sym"],c=SC[r["scen"]][1],marker=VAR[r["var"]][1],s=40,edgecolors="k",lw=0.4)
ax.set_xscale("log"); ax.axvline(1,color="green",lw=1.5); ax.axhline(0.5,color="0.6",ls=":",lw=0.8); ax.axhline(0.9,color="0.4",ls="--",lw=0.8)
ax.text(1.05,0.91,"near-symmetric (≥90%)",fontsize=7,color="0.3"); ax.text(1.05,0.51,"no preference (50%)",fontsize=7,color="0.4"); ax.text(1.05,0.02,"1× CMB level",color="green",fontsize=7)
ax.set_xlim(0.7,5000); ax.set_ylim(0,1.02); ax.set_xlabel("final leftover (× floor, log scale)"); ax.set_ylabel("share of residual power that is mirror-symmetric\nabout the CoE axis (ℓ 15°/195°), modes 1–12")
ax.set_title("Mirror symmetry of the final residual, tidal-depth fill, core empty (one point per seed)")
hd=[Line2D([],[],marker="o",ls="",color=c,ms=7,label=n) for n,c in SC.values()]+[Line2D([],[],marker=mk,ls="",mfc="white",mec="k",ms=7,label=n) for n,mk in VAR.values()]
ax.legend(handles=hd,loc="upper center",bbox_to_anchor=(0.5,-0.12),ncol=4,frameon=False,fontsize=7.5)
fig.tight_layout(); fig.savefig(f"{OUT}/symmetry_tidal_vs_leftover.png",dpi=200,bbox_inches="tight"); plt.close(fig)
# chart 2: per-mode symmetric share, median over seeds, one panel per variant
fig,axs=plt.subplots(1,1,figsize=(7.1,3.6)); axs=[axs]
for ax,(v,(vn,_)) in zip(axs,VAR.items()):
    for k,(sc,(n,c)) in enumerate(SC.items()):
        S=[r["sym_m"] for r in rows if r["scen"]==sc and r["var"]==v]
        if not S: continue
        S=np.array(S); ax.bar(m+(k-1.5)*0.2,np.median(S,0),width=0.2,color=c,label=f"{n} ({len(S)} seeds)")
        ax.vlines(m+(k-1.5)*0.2,S.min(0),S.max(0),color="k",lw=0.6)
    ax.axhline(0.5,color="0.5",ls=":",lw=0.8); ax.set_ylim(0,1); ax.set_ylabel("symmetric share"); ax.set_title(vn)
axs[-1].set_xticks(m); axs[-1].set_xlabel("mode m")
fig.legend(handles=[plt.Rectangle((0,0),1,1,fc=c,label=n) for n,c in SC.values()],loc="lower center",ncol=4,frameon=False,fontsize=7.5)
fig.text(0.5,0.035,"bars = median over seeds, black lines = range across seeds; 0.5 = no preference",ha="center",fontsize=7)
fig.tight_layout(rect=(0,0.06,1,1)); fig.savefig(f"{OUT}/symmetry_tidal_by_mode.png",dpi=200,bbox_inches="tight"); plt.close(fig)
for sc in SC:
    for v in VAR:
        L=[r for r in rows if r["scen"]==sc and r["var"]==v]
        if L: print(sc,v,len(L),"left",[round(r["left"]) for r in sorted(L,key=lambda r:r["seed"])],"sym",[round(r["sym"],2) for r in sorted(L,key=lambda r:r["seed"])],"body_sym",[round(r["body_sym"],2) for r in L],"sym rms",[round(r["sym_rms"]) for r in L],"anti rms",[round(r["antisym_rms"]) for r in L])

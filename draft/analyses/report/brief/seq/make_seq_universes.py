"""Final universes of the six sequential-placement curves of Fig. 7 (side view + leftover around the ring).
A: saved layouts (../../../sequential/local_frozen.npz, shallow/near_fit.json). B: regenerated with regen_seq.py (same code+seed)."""
import os,sys,json,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
HERE=os.path.dirname(os.path.abspath(__file__)); FIG=os.path.join(HERE,"..","figures"); SEQ="/workspace/clump_scan/sequential"
sys.path.insert(0,SEQ); _c=os.getcwd(); os.chdir(SEQ); from seq_core import R,FLOOR,rms,C_G,A_IN_G,RMIN,S_MAX_G; os.chdir(_c)
plt.rcParams.update({"font.size":7.5,"axes.titlesize":7.5,"savefig.dpi":300,"pdf.fonttype":42})
STY={"seqA_deep":"#c0392b","seqA_shallow":"#e59866","seqB_cmb_aware":"#7d3c98","seqB_foam_rim":"#5dade2","seqB_random":"#7f8c8d","seqB_shallow":"#d7bde2"}
NAME={"seqA_deep":"A: near fix, deep voids + clusters","seqA_shallow":"A: near fix, shallow voids","seqB_cmb_aware":"B: growth, CMB-aware",
      "seqB_foam_rim":"B: growth, foam rule","seqB_random":"B: growth, random","seqB_shallow":"B: growth, CMB-aware, shallow voids"}
SLAB=1.5   # B runs: bodies whose centre lies within 1.5 Gpc of the ring plane (their 3D layout fills the whole shell region)
E={}
z=np.load(f"{SEQ}/local_frozen.npz",allow_pickle=True); k=z["kind"].astype(str)
E["seqA_deep"]=dict(Xv=z["X"][k=="void"],rv=z["r"][k=="void"],Xc=z["X"][k=="clu"],Mc=z["M"][k=="clu"],w=z["w"],n=f"{(k=='void').sum()} deep voids + {(k=='clu').sum()} clusters",slab=False)
b=json.load(open(f"{SEQ}/shallow/near_fit.json"))["best"]
E["seqA_shallow"]=dict(Xv=np.array(b["X"]),rv=np.array(b["Rv"]),Xc=np.zeros((0,3)),Mc=np.zeros(0),w=np.array(b["w"]),n=f"{b['n']} shallow voids",slab=False)
for key,rule,seed in [("seqB_cmb_aware","cmb_aware",1),("seqB_foam_rim","foam_rim",3),("seqB_random","random",2)]:
    fn=f"{HERE}/final_b_{rule}_s{seed}.npz"
    if not os.path.exists(fn) and os.environ.get("PREVIEW"): fn=f"{HERE}/final_b_foam_rim_s3.npz"
    d=np.load(fn); X=d["Xp"]; r=d["rp"]; kp=d["kp"]; v=kp=="void"
    E[key]=dict(Xv=X[v],rv=r[v],Xc=X[~v],Mc=np.ones((~v).sum()),w=d["w"],n=f"{v.sum():,} voids + {(~v).sum():,} clusters",slab=True,chk=(float(d["left"]),float(d["ref"])))
d=np.load(f"{HERE}/final_shallow_cmb_aware_s2.npz"); X=d["X"]
E["seqB_shallow"]=dict(Xv=X,rv=np.full(len(X),1.5),Xc=np.zeros((0,3)),Mc=np.zeros(0),w=d["w"],n=f"{len(X):,} shallow voids",slab=True,chk=(float(d["left"]),float(d["ref"])))
for kk,e in E.items(): e["left"]=float(rms(R-e["w"])/FLOOR)
order=sorted(E,key=lambda q:-E[q]["left"])
for q in order: print(q,round(E[q]["left"],2),E[q].get("chk"),E[q]["n"])
json.dump({q:dict(left=E[q]["left"],n=E[q]["n"],check=E[q].get("chk")) for q in order},open(f"{HERE}/seq_universes.json","w"),indent=1)
# Chris: same format as the example-universe figures (one row per run: side view | leftover by mode | residual), 3 rows per figure
plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"axes.labelsize":8,"xtick.labelsize":7,"ytick.labelsize":7})
TH=np.radians(np.arange(360))
def amps(v,M=12):
    m=np.arange(1,M+1)[:,None]; return np.hypot(2*np.mean(v[None]*np.cos(m*TH),1),2*np.mean(v[None]*np.sin(m*TH),1))
def hi13(v):
    rec=np.full(360,v.mean())
    for m in range(1,13): rec+=2*np.mean(v*np.cos(m*TH))*np.cos(m*TH)+2*np.mean(v*np.sin(m*TH))*np.sin(m*TH)
    return float(np.sqrt(2)*rms(v-rec)/FLOOR)
START=amps(R)/FLOOR; W=46; Ld=np.arange(360)
def row(axs,i,q):
    e=E[q]; col=STY[q]; tcol="#9b6fb5" if q=="seqB_shallow" else col   # pale curve colour darkened for readable text
    ax=axs[0]
    ax.add_patch(Circle(C_G[:2],A_IN_G,fill=False,ec="0.6",lw=0.8,ls="--",zorder=1))
    Xv,rv,Xc=e["Xv"],e["rv"],e["Xc"]
    if e["slab"]: sv=np.abs(Xv[:,2])<SLAB; sc=np.abs(Xc[:,2])<SLAB; Xv,rv,Xc=Xv[sv],rv[sv],Xc[sc]
    ax.add_collection(PatchCollection([Circle(x[:2],r) for x,r in zip(Xv,rv)],fc="#a9cce3",ec="#5d8aa8",lw=0.25,alpha=0.9,zorder=2))
    if len(Xc): ax.scatter(Xc[:,0],Xc[:,1],s=1.0 if e["slab"] else 7,c="#c0392b",edgecolors="k" if not e["slab"] else "none",linewidths=0.3,zorder=3)
    ax.add_patch(Circle((0,0),RMIN,fill=False,ec="k",lw=1.3,zorder=4)); ax.plot(0,0,"k+",ms=8,mew=1.2,zorder=5)
    ax.set_xlim(-W,W); ax.set_ylim(-W,W); ax.set_aspect("equal"); ax.set_xticks([-40,-20,0,20,40]); ax.set_yticks([-40,-20,0,20,40])
    ax.tick_params(labelsize=6); ax.set_xlabel("Gpc",fontsize=7); ax.set_ylabel("Gpc",fontsize=7)
    ax.set_title(f"{i+1}. {NAME[q].replace("growth, CMB-aware, shallow voids","CMB-aware, shallow voids")}: {e['left']:,.1f}×\n{e['n']}",fontsize=7.5,color=tcol)
    v=R-e["w"]; a=axs[1]; mm=np.arange(1,13); a.bar(mm,np.maximum(amps(v)/FLOOR,1e-4),color=col)
    a.bar(mm,START,color="none",edgecolor="0.35",lw=0.6,ls="--",zorder=3); a.bar(13.8,max(hi13(v),1e-4),color="0.55",width=1.3)
    a.axhline(1,color="green",lw=1.2); a.set_yscale("log"); a.set_ylim(1e-3,4000); a.set_xticks([1,3,5,7,9,11,13.8]); a.set_xticklabels(["1","3","5","7","9","11","13+"])
    a.set_xlabel("mode $m$"); a.set_ylabel(r"$\times$ floor"); a.set_title("leftover by mode (dashed = start)",fontsize=7.5)
    c=axs[2]; res=v/FLOOR; c.axhline(0,color="0.6",lw=0.5,zorder=1)
    c.axhspan(-1,1,fc=(0.0,0.5,0.0,0.25),ec=(0.0,0.5,0.0,0.9),lw=0.6,zorder=1.5); c.plot(Ld,res,color=col,lw=0.9,zorder=2)
    lim=max(1.15*np.abs(res).max(),1.05); c.set_ylim(-lim,lim)
    c.set_xlim(0,360); c.set_xticks(range(0,361,90)); c.set_xlabel(r"$\ell$ (deg)"); c.set_ylabel(r"$\times10^{-5}$")
    c.set_title(f"residual after placement: {e['left']:,.1f}×",fontsize=7.5,color=tcol)
for part,rows in (("a",order[:3]),("b",order[3:])):
    fig,axs=plt.subplots(3,3,figsize=(7.1,2.3*3),gridspec_kw=dict(width_ratios=[0.9,1,1.25]))
    for j,q in enumerate(rows): row(axs[j],order.index(q),q)
    fig.tight_layout(h_pad=1.0,w_pad=1.0)
    fig.savefig(f"{FIG}/fig_seq_universes_{part}.pdf",bbox_inches="tight"); fig.savefig(f"{FIG}/fig_seq_universes_{part}.png",bbox_inches="tight",dpi=200); plt.close(fig)

"""Triangle area histograms, full Delaunay (no cutoff), foam s2: BEFORE grey vs AFTER dark blue, log bins; whole slab and 0-3 Gpc band (by centroid)."""
import json,numpy as np,matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
src=open("tri_view.py").read(); exec(src[:src.index("def links(")])
d=json.load(open("minedit_ent/b_foam_rim_s2_tidalcore_entropy.json")); base=d["tag"].replace("_entropy","")
out={}
for k,L in (("before",EB.layout(base)),("after",EB.layout(base,d))):
    sl,P,S,E,A=tri(L[0],L[1]); c=P[S].mean(1); out[k]=(A,np.abs(np.hypot(c[:,0],c[:,1])-RMIN))
bins=np.logspace(-3,np.log10(max(out["before"][0].max(),out["after"][0].max())*1.05),36); info={}
for tag,sel,fn in (("whole slab",lambda dr: np.ones(len(dr),bool),"tri_foam_s2_area_hist.png"),("0–3 Gpc from the CMB ring",lambda dr: dr<3,"tri_foam_s2_area_hist_ring.png")):
    fig,ax=plt.subplots(figsize=(8.5,5.6)); info[tag]={}
    for k,c,ls in (("before","0.55","-"),("after","#1f3a93","-")):
        A=out[k][0][sel(out[k][1])]; info[tag][k]=dict(n=int(len(A)),median=float(np.median(A)))
        ax.hist(A,bins,histtype="step",color=c,lw=2.2,label=f"{k.upper()}: {len(A):,} triangles, median {np.median(A):.2f} Gpc²")
    ax.set_xscale("log"); ax.set_xlabel("triangle area (Gpc², log bins)",fontsize=12); ax.set_ylabel("number of triangles",fontsize=12); ax.tick_params(labelsize=11); ax.grid(alpha=0.3)
    ax.set_title(f"Triangle areas, full Delaunay triangulation — foam s2, {tag}",fontsize=12.5)
    ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.15),ncol=1,frameon=False,fontsize=11.5)
    fig.savefig(fn,dpi=140,bbox_inches="tight"); plt.close(fig)
print(json.dumps(info,indent=1))

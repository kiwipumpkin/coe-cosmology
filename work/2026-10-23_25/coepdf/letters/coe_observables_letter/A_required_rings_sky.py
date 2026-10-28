# Sky maps of the best free-axis and best shared-axis ring solutions, ring axes marked. From A_required_rings_results.json.
import numpy as np,matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
from A_required_rings_common import *
best=min([k for k in R if k.startswith("free")],key=lambda k:R[k]["resid"]);sh=min([k for k in R if k.startswith("shared_fit")],key=lambda k:R[k]["resid"])
fig,ax=plt.subplots(1,2,figsize=(15,5))
for a,k in zip(ax,[sh,best]):
  s=sky(R[k]["rings"]).reshape(LL.shape);v=np.abs(s).max();im=a.imshow(s,origin="lower",extent=[0,360,-90,90],cmap="RdBu_r",vmin=-v,vmax=v,aspect="auto");plt.colorbar(im,ax=a,label="delta (+ matter, - void)")
  for r in R[k]["rings"]:a.plot(r[0]%360,r[1],"k*" ,ms=10,mfc="yellow")
  a.plot(180,-15,"ms",label="CoE");a.plot(260,60,"g^",label="AoE");a.set(xlabel="l [deg]",ylabel="b [deg]",title=f"{k}: residual rms {R[k]['resid']:.4f}");a.legend(fontsize=8)
fig.text(.5,.01,"Yellow stars: ring axes. Not unique: latitude structure unconstrained by the strip average.",ha="center",fontsize=9)
plt.tight_layout(rect=(0,.04,1,1));plt.savefig("A_required_rings_sky.png",dpi=140)

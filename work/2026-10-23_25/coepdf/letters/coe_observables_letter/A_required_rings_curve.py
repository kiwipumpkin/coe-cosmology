# A/A_all before/after adding strip-averaged ring delta, with Migkas data. From A_required_rings_results.json.
import numpy as np,matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
from A_required_rings_common import *
best=min([k for k in R if k.startswith("free")],key=lambda k:R[k]["resid"]);sh=min([k for k in R if k.startswith("shared_fit")],key=lambda k:R[k]["resid"])
plt.figure(figsize=(11,5));plt.errorbar(l,y,yerr=[elo,ehi],fmt="k.",ms=4,lw=.7,label="Migkas Fig.6");plt.plot(l,y,"k-",lw=.8,label="before")
for k,c in [(sh,"b"),(best,"r")]:D=St@sky(R[k]["rings"]);plt.plot(l,y+D-D.mean(),c+".-",lw=1,label=f"after {k} (mean removed)")
plt.xlabel("l [deg]");plt.ylabel("A/A_all");plt.xlim(0,360);plt.grid(alpha=.3);plt.legend(fontsize=8);plt.tight_layout();plt.savefig("A_required_rings_curve.png",dpi=140)

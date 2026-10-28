# Residual rms (m=2..35) and delta rms vs number of rings for each axis case; blob solution as reference line. From A_required_rings_results.json.
import numpy as np,matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
from A_required_rings_common import *
best=min([k for k in R if k.startswith("free")],key=lambda k:R[k]["resid"]);sh=min([k for k in R if k.startswith("shared_fit")],key=lambda k:R[k]["resid"])
fig,ax=plt.subplots(1,2,figsize=(13,4.5))
for c,col in zip(["CoE","AoE","shared_fit","free"],"mgbr"):
  N=[n for n in [2,4,6,8] if f"{c}_{n}" in R];ax[0].semilogy(N,[R[f"{c}_{n}"]["resid"] for n in N],col+"o-",label=c);ax[1].semilogy(N,[R[f"{c}_{n}"]["drms"] for n in N],col+"o-",label=c)
ax[0].axhline(0.1272,color="k",ls=":",label="no change (0.127)");ax[0].axhline(0.0058,color="k",ls="--",label="blob/delta(l) solution (0.0058)");ax[1].axhline(1.41,color="k",ls="--",label="delta(l) solution rms 1.41")
ax[0].set(xlabel="number of rings",ylabel="residual rms m=2..35");ax[1].set(xlabel="number of rings",ylabel="delta rms");[a.legend(fontsize=8) or a.grid(alpha=.3,which="both") for a in ax]
plt.tight_layout();plt.savefig("A_required_rings_vsN.png",dpi=140)

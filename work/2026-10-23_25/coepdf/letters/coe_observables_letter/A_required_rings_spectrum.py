# Shares m=2..35 (dipole removed) vs Planck and residuals for best shared-axis and free-axis ring fits. From A_required_rings_results.json.
import numpy as np,matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
from A_required_rings_common import *
best=min([k for k in R if k.startswith("free")],key=lambda k:R[k]["resid"]);sh=min([k for k in R if k.startswith("shared_fit")],key=lambda k:R[k]["resid"])
fig,ax=plt.subplots(1,2,figsize=(14,5));a0=share(y)
ax[0].semilogy(m,pl,"rs-",ms=4,label="Planck");ax[0].semilogy(m,a0,"k.-",label="before")
ax[1].plot(m,a0-pl,"k.-",label=f"before rms {np.sqrt(((a0-pl)**2).mean()):.4f}")
for k,c in [(sh,"b"),(best,"g")]:a=share(y+St@sky(R[k]["rings"]));ax[0].semilogy(m,a,c+".-",label=k);ax[1].plot(m,a-pl,c+".-",label=f"{k} rms {np.sqrt(((a-pl)**2).mean()):.4f}")
ax[0].set(xlabel="m = l",ylabel="share",title="Shares m=2..35");ax[0].grid(alpha=.3,which="both");ax[0].legend(fontsize=8)
ax[1].axhline(0,color="k",lw=.6);[ax[1].axvline(q,color="gray",lw=.5,ls=":") for q in range(4,36,4)];ax[1].set(xlabel="m = l",ylabel="Migkas - Planck share",title="Residuals (dotted m=4k)");ax[1].grid(alpha=.3);ax[1].legend(fontsize=8)
fig.text(.5,.01,"Dipole removed. m azimuthal, l spherical: not the same quantity.",ha="center",fontsize=9);plt.tight_layout(rect=(0,.04,1,1));plt.savefig("A_required_rings_spectrum.png",dpi=140)

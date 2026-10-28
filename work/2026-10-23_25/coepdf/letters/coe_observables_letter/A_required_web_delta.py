# (a) Required delta: latitude-independent minimal-norm delta(l) and the cluster/void NNLS sky map. Data from A_required_web_solve.py.
import numpy as np,matplotlib;matplotlib.use("Agg");import matplotlib.pyplot as plt
R=np.load("A_required_web_results.npz",allow_pickle=True)
fig,ax=plt.subplots(1,2,figsize=(15,5))
ax[0].plot(R["lg"],R["x"],"k-",lw=1);ax[0].axhline(0,color="gray",lw=.6);ax[0].fill_between(R["lg"],R["x"],0,where=R["x"]>0,color="r",alpha=.3,label="matter (+)");ax[0].fill_between(R["lg"],R["x"],0,where=R["x"]<0,color="b",alpha=.3,label="voids (-)")
ax[0].set(xlabel="l [deg]",ylabel="delta (fractional A/A_all units)",title=f"Latitude-independent delta(l): rms {R['x'].std():.2f}, peak {abs(R['x']).max():.2f}",xlim=(0,360));ax[0].legend();ax[0].grid(alpha=.3)
mp=R["mapb"];v=np.percentile(abs(mp),99);im=ax[1].imshow(mp,origin="lower",extent=[0,360,-90,90],cmap="RdBu_r",vmin=-v,vmax=v,aspect="auto");plt.colorbar(im,ax=ax[1],label="delta")
km=R["keep_meta"];kc=km[:,0]=="cluster"
ax[1].scatter(km[kc,2].astype(float),km[kc,3].astype(float),s=8,c="darkred",marker="+",label=f"clusters 3-8 deg: {kc.sum()}")
ax[1].scatter(km[~kc,2].astype(float),km[~kc,3].astype(float),s=10,facecolors="none",edgecolors="navy",label=f"voids 15-30 deg: {(~kc).sum()}")
ax[1].set(xlabel="l [deg]",ylabel="b [deg]",title="Cluster/void version (NNLS, candidates on 10 deg x 30 deg grid)");ax[1].legend(fontsize=8,loc="lower right")
fig.text(.5,.01,"Not unique: phases free and latitude structure unconstrained by the 90 deg strip average. Colour scale clipped at 99th percentile.",ha="center",fontsize=9)
plt.tight_layout(rect=(0,.04,1,1));plt.savefig("A_required_web_delta.png",dpi=140)

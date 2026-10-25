import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
ri=np.arange(3.0,9.0); ro=ri+1; m=(ri+ro)/2; V=(4/3)*np.pi*(ro**3-ri**3)   # innermost first
p=brentq(lambda q:(m[0]**q/V[0])/(m[-1]**q/V[-1])-1.03,0.5,4)
assert abs(p-1.9603)<5e-5
S={"(a) equal mass per shell":np.ones(6),"(b) mass ∝ volume (flat density)":V/V[0],f"(c) mass ∝ r^{p:.4f} (3% density drop)":(m/m[0])**p}
rho={k:q/V for k,q in S.items()}
kb,kc=list(S)[1],list(S)[2]
assert np.ptp(rho[kb])<1e-15 and abs(rho[kc][0]/rho[kc][-1]-1.03)<1e-12 and np.all(np.diff(rho[kc])<0)
fig,ax=plt.subplots(1,2,figsize=(16,6.5),sharex=True)
for (k,q),mk,col in zip(S.items(),"osD",["tab:red","tab:gray","tab:purple"]):
    edges=np.r_[ri,ro[-1]]
    ax[0].stairs(q,edges,color=col,lw=2.5,label=k,baseline=None); ax[1].stairs(rho[k]/rho[k][0],edges,color=col,lw=2.5,label=k,baseline=None)
ax[0].set_ylabel("mass (matter + light) per shell, innermost shell = 1")
ax[1].set_ylabel("density per shell ÷ innermost shell density")
ax[0].set_title("Mass per shell (constant within each shell)",pad=18); ax[1].set_title("Density per shell (mass ÷ shell volume)",pad=18)
for a in ax:
    a.axvline(3,color="k",lw=2); a.axvline(9,color="k",lw=2); a.set_xlim(3,9); a.grid(alpha=.3)
    a.set_xticks(range(3,10)); a.set_xticklabels(["D−R\n3","4","5","6","7","8","D+R\n9"])
    for e in range(4,9): a.axvline(e,color="0.85",lw=0.8)
    for k_,mm_ in zip(range(6,0,-1),m): a.text(mm_,1.01,f"shell {k_}",transform=a.get_xaxis_transform(),ha="center",fontsize=9)
    a.set_xlabel("distance from CoE (c·Δt); shell edges at integers")
fig.legend(*ax[0].get_legend_handles_labels(),loc="lower center",ncol=3,fontsize=11,frameon=False)
fig.suptitle("Exactly-c geometry, t = 3: D = 6, R = 3, six shells of width c·Δt (spherical volumes)",fontsize=13)
plt.tight_layout(rect=[0,0.08,1,0.95]); plt.savefig("scenarios_mass_density.png",dpi=180)
print("p",p); print("shell r_mid V", *[f"{x:.3f}" for x in V])
for k in S: print(k,"mass",np.round(S[k],4),"density",np.round(rho[k]/rho[k][0],4))

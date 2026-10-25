exec(open("paper_repro/stretch/gr_exact.py").read().split("p1=np.array")[0])
import json
fullN,fullG=np.load("paper_repro/stretch/gr_top_curves.npy")
print("csv vs recomputed Newtonian max diff",np.abs(fullN-full).max())
st=json.load(open("/tmp/gr_top_stats.json"))
GN=fullN/fullN.mean(); GG=fullG/fullG.mean()
p1=np.array(P2full,float); p0=p1.copy(); p0[2]=0
out={}
for nm,p,fn,G in (("weak-field, mass (Newtonian geometry)",p1,S_model,GN),("exact GR, mass",p1,S_gr,GG),("exact GR, zero mass",p0,S_gr,GG)):
    o=fn(*p); out[nm]=resid(G,o[0]); print(nm,rms(out[nm]))
b_=np.loadtxt("migkas-A-vs-angle-bins.csv",delimiter=",",skiprows=1)
plt.rcParams.update({"font.size":18}); fig,ax=plt.subplots(2,1,figsize=(13,12),sharex=True)
ax[0].errorbar(b_[:,0],b_[:,1],yerr=[b_[:,2],b_[:,3]],fmt="o",color="k",ms=5,capsize=2,label="Migkas bins")
ax[0].plot(ELL,fullN,lw=3,color="tab:orange",label=f"Newtonian: χ² {st['chiN'][0]:.1f}")
ax[0].plot(ELL,fullG,lw=3,ls="--",color="tab:blue",label=f"exact GR: χ² {st['chiG'][0]:.1f}")
ax[0].set_ylabel("A / A_all"); ax[0].legend(fontsize=14,loc="upper center",bbox_to_anchor=(0.5,1.18),ncol=3); ax[0].grid(alpha=.3)
for (nm,R),c,ls in zip(out.items(),("tab:red","tab:blue","tab:green"),("-","--","-")):
    ax[1].plot(ELL,R*1e3,color=c,ls=ls,lw=3,label=f"{nm}: RMS {rms(R):.2e}")
ax[1].axhline(0,color="gray"); ax[1].set_ylabel("CMB ΔT/T × 10³"); ax[1].set_xlabel("galactic longitude ℓ (deg)")
ax[1].legend(fontsize=14,loc="upper center",bbox_to_anchor=(0.5,-0.2),ncol=1); ax[1].grid(alpha=.3); ax[1].set_xlim(0,359)
plt.tight_layout(); plt.savefig("paper_repro/stretch/gr_exact_both.png",dpi=120,bbox_inches="tight")

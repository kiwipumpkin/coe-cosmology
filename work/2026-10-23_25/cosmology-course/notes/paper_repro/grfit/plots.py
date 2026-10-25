from core import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
t0=time.time()
def best(f): d=json.load(open(f)); return min(d["rows"],key=lambda r:r["chi2"])
g=best("fitb_ltb_plain.json"); gs=best("fitb_ltb_sine.json"); n=best("fit_newt_plain.json")
cg=shifted(ltb(g["D"],g["mcc"],g["kin"],g["kout"]),g["l0"])
cs=sine(shifted(ltb(gs["D"],gs["mcc"],gs["kin"],gs["kout"]),gs["l0"]),gs["l0"])[0]
cn=shifted(nsky(n["D"],n["mcc"],n["kin"],n["kout"]),n["l0"])
print("check",chi(cg),chi(cs),chi(cn))
b_=np.loadtxt("../../migkas-A-vs-angle-bins.csv",delimiter=",",skiprows=1)
plt.rcParams.update({"font.size":18}); fig,ax=plt.subplots(figsize=(13,8))
ax.errorbar(b_[:,0],b_[:,1],yerr=[b_[:,2],b_[:,3]],fmt="o",color="k",ms=5,capsize=2,label="Migkas bins")
ax.plot(TH1,cn,lw=3,color="tab:orange",label=f"matched-today GR (= Newtonian): χ² {chi(cn):.1f}")
ax.plot(TH1,cg,lw=3,ls="--",color="tab:blue",label=f"nested-shell GR: χ² {chi(cg):.1f}")
ax.plot(TH1,cs,lw=3,color="tab:green",label=f"nested-shell GR + sine: χ² {chi(cs):.1f}")
ax.set_xlim(0,359); ax.set_xlabel("galactic longitude ℓ (deg)"); ax.set_ylabel("A / A_all"); ax.grid(alpha=.3)
ax.legend(fontsize=15,loc="lower center",bbox_to_anchor=(0.5,1.01),ncol=2)
plt.tight_layout(); plt.savefig("../stretch/gr_bestfit.png",dpi=120,bbox_inches="tight")
# profile D vs mcc, nested GR no sine, kin,kout fixed at best, l0 minimised on 1-deg grid
Ds=np.geomspace(20,150,36); Ms=np.geomspace(1.7e-4,0.17,36); C=np.full((36,36),np.nan)
for i,D in enumerate(Ds):
    for j,m in enumerate(Ms):
        af=ltb(D,m,g["kin"],g["kout"])
        if af is None: continue
        C[j,i]=min(chi(np.roll(af,k)) for k in range(360))
np.save("profile.npy",C)
fig,ax=plt.subplots(figsize=(13,8)); dc=C-np.nanmin(C)
pc=ax.pcolormesh(Ds,Ms,np.minimum(dc,25),shading="auto",cmap="viridis"); ax.contour(Ds,Ms,dc,[1,4],colors="w")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("D (Gpc)"); ax.set_ylabel("mcc"); fig.colorbar(pc,label="Δχ² (capped at 25)")
ax.set_title("nested-shell GR, no sine; white lines Δχ²=1, 4")
plt.tight_layout(); plt.savefig("../stretch/gr_bestfit_degeneracy.png",dpi=120)
k=np.unravel_index(np.nanargmin(C),C.shape); print("profile min",np.nanmin(C),"D",Ds[k[1]],"mcc",Ms[k[0]],"cells within 1:",int((dc<1).sum()),"Drange",Ds[(dc<1).any(0)].min(),Ds[(dc<1).any(0)].max(),"mrange",Ms[(dc<1).any(1)].min(),Ms[(dc<1).any(1)].max(),"t",time.time()-t0)

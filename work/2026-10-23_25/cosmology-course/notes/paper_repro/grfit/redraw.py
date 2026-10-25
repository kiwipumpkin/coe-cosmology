from core import *
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
def pick(f,D):
    return min([r for r in json.load(open(f))["rows"] if abs(r["D"]-D)<0.1],key=lambda r:r["chi2"])
n=pick("fit_newt_plain.json",58.1967); g=pick("fit_ltb_plain.json",55.4648); s=pick("fit_ltb_sine.json",42.3144)
cn=shifted(nsky(n["D"],n["mcc"],n["kin"],n["kout"]),n["l0"])
cg=shifted(ltb(g["D"],g["mcc"],g["kin"],g["kout"]),g["l0"])
cs=sine(shifted(ltb(s["D"],s["mcc"],s["kin"],s["kout"]),s["l0"]),s["l0"])[0]
for nm,r,c in (("newt",n,cn),("gr",g,cg),("grs",s,cs)): print(nm,r["D"],r["mcc"],r["l0"],"saved",r["chi2"],"recomputed",chi(c))
b_=np.loadtxt("../../migkas-A-vs-angle-bins.csv",delimiter=",",skiprows=1)
plt.rcParams.update({"font.size":18}); fig,ax=plt.subplots(figsize=(13,8))
ax.errorbar(b_[:,0],b_[:,1],yerr=[b_[:,2],b_[:,3]],fmt="o",color="k",ms=5,capsize=2,label="Migkas bins")
ax.plot(TH1,cn,lw=3,color="tab:orange",label=f"Newtonian, D {n['D']:.1f} Gpc: χ² {chi(cn):.2f}")
ax.plot(TH1,cg,lw=3,ls="--",color="tab:blue",label=f"nested-shell GR, D {g['D']:.1f} Gpc: χ² {chi(cg):.2f}")
ax.plot(TH1,cs,lw=3,color="tab:green",label=f"nested-shell GR with sine, D {s['D']:.1f} Gpc: χ² {chi(cs):.2f}")
ax.set_xlim(0,359); ax.set_xlabel("galactic longitude ℓ (deg)"); ax.set_ylabel("A / A_all"); ax.grid(alpha=.3)
ax.legend(fontsize=15,loc="lower center",bbox_to_anchor=(0.5,1.01),ncol=2)
plt.tight_layout(); plt.savefig("../stretch/gr_bestfit.png",dpi=120,bbox_inches="tight")

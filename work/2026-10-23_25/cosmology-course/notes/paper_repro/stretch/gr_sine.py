exec(open("paper_repro/stretch/gr_exact.py").read().split("p1=np.array")[0])
import sys,json; sys.path.insert(0,"/workspace/cosmology/papers/spreadsheet_model/integrated")
from integ_core import best_sine, A_from
from integ_core import ML, MA, MLO, MHI; IMIG=ML.astype(int)
fullN,fullG=np.load("paper_repro/stretch/gr_top_curves.npy"); l0=15.1629523788697
th=ELL-l0; amp,phi=best_sine(fullG[IMIG],th[IMIG],MA)
fullS=A_from(fullG,th,amp,phi)
W=1/(0.5*(MLO+MHI))**2
def st(m): mm=m[IMIG]; return float(np.sum(W*(mm-MA)**2)),float(1-np.sum((MA-mm)**2)/np.sum((MA-MA.mean())**2))
s0,s1=st(fullG),st(fullS); print("amp",amp,"phi deg",np.degrees(phi),"before",s0,"after",s1)
GN=fullS/fullS.mean(); print("RMS(G-1) sine",rms(GN-1))
p1=np.array(P2full,float); p0=p1.copy(); p0[2]=0; out={}
for nm,p,fn in (("weak-field, mass",p1,S_model),("exact GR, mass",p1,S_gr),("exact GR, zero mass",p0,S_gr)):
    o=fn(*p); out[nm]=resid(GN,o[0]); print(nm,rms(out[nm]))
json.dump(dict(amp=amp,phi_deg=float(np.degrees(phi)),before=s0,after=s1,rmsG=rms(GN-1),rms={k:rms(v) for k,v in out.items()}),open("paper_repro/stretch/gr_sine.json","w"),indent=1)
b_=np.loadtxt("migkas-A-vs-angle-bins.csv",delimiter=",",skiprows=1)
plt.rcParams.update({"font.size":18}); fig,ax=plt.subplots(2,1,figsize=(13,12),sharex=True)
ax[0].errorbar(b_[:,0],b_[:,1],yerr=[b_[:,2],b_[:,3]],fmt="o",color="k",ms=5,capsize=2,label="Migkas bins")
ax[0].plot(ELL,fullG,lw=3,color="tab:orange",label=f"no sine: χ² {s0[0]:.1f}")
ax[0].plot(ELL,fullS,lw=3,ls="--",color="tab:blue",label=f"with sine: χ² {s1[0]:.1f}")
ax[0].set_ylabel("A / A_all"); ax[0].legend(fontsize=14,loc="lower center",bbox_to_anchor=(0.5,1.01),ncol=3); ax[0].grid(alpha=.3)
for (nm,R),c,ls in zip(out.items(),("tab:red","tab:blue","tab:green"),("-","--","-")):
    ax[1].plot(ELL,R*1e3,color=c,ls=ls,lw=3,label=f"{nm}: RMS {rms(R):.2e}")
ax[1].axhline(0,color="gray"); ax[1].set_ylabel("CMB ΔT/T × 10³"); ax[1].set_xlabel("galactic longitude ℓ (deg)")
ax[1].legend(fontsize=14,loc="upper center",bbox_to_anchor=(0.5,-0.2)); ax[1].grid(alpha=.3); ax[1].set_xlim(0,359)
plt.tight_layout(); plt.savefig("paper_repro/stretch/gr_sine.png",dpi=120,bbox_inches="tight")

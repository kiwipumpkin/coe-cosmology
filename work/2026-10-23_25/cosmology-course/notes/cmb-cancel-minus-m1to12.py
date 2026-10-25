import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src=open("cmb-cancel.py").read().split('if __name__')[0]; g={}; exec(src,g)
out=json.load(open("cmb-cancel.json")); ELL=g["ELL"]; th=np.radians(ELL)
def strip(R,M=12):
    R=R-R.mean()
    for m in range(1,M+1):
        a=2*np.mean(R*np.cos(m*th)); b=2*np.mean(R*np.sin(m*th)); R=R-a*np.cos(m*th)-b*np.sin(m*th)
    return R
def P(k): r=out[k]; return [r[x] for x in ("beta1","q","mu","s","phi_deg")]
cases=[("run2","s ≤ 3, full precision",g["G_FULL"],"C0"),("run3","s ≤ 50, full precision",g["G_FULL"],"C2"),
       ("run2_rounded","s ≤ 3, 4-decimal CSV",g["G_RND"],"C3")]
plt.rcParams.update({"font.size":15})
fig,ax=plt.subplots(2,1,figsize=(12,10),sharex=True)
res={}
for k,lab,G,c in cases:
    L=strip(g["resid"](G,g["S_model"](*P(k))[0])); rm=float(np.sqrt(np.mean(L**2))); mx=float(np.abs(L).max())
    res[k]=dict(rms=rm,max=mx)
    a=ax[1] if "rounded" in k else ax[0]
    a.plot(ELL,L*1e5,color=c,lw=1.8 if "rounded" not in k else 1.2,label=f"{lab}: RMS {rm:.2e}, max {mx:.2e}")
for a,t in zip(ax,["Full-precision curve","Rounded 4-decimal curve (as in your plot)"]):
    a.axhline(1,ls="--",color="gray");a.axhline(-1,ls="--",color="gray");a.axhline(0,color="k",lw=.6)
    a.set_ylabel(r"Residual $\times10^{5}$");a.grid(alpha=.3);a.set_title(t,pad=75 if a is ax[0] else 48)
    a.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=1,frameon=False,fontsize=14)
ax[1].set_xlabel(r"Galactic longitude $\ell$ (deg)");ax[1].set_xlim(0,360)
fig.suptitle(r"Residual after removing $m=1$ through $12$ (dashed lines at $\pm10^{-5}$)",fontsize=17)
fig.tight_layout(); fig.savefig("cmb-cancel-minus-m1to12.png",dpi=130)
json.dump(res,open("cmb-cancel-minus-m1to12.json","w"),indent=1); print(res)

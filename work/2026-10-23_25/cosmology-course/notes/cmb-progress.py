import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src=open("cmb-cancel.py").read().split('if __name__')[0]; g={}; exec(src,g)
d=json.load(open("cmb-cancel-story.json")); P=d["params"]; phi=P[4]
ELL=g["ELL"]; G=g["G_FULL"]; th=np.radians(ELL-phi); S=g["S_model"](*P)[0]; R=g["resid"](G,S)
def strip(x,M):
    x=x-x.mean()
    for m in range(1,M+1):
        a=2*np.mean(x*np.cos(m*th)); b=2*np.mean(x*np.sin(m*th)); x=x-a*np.cos(m*th)-b*np.sin(m*th)
    return x
rms=lambda x: float(np.sqrt(np.mean(x**2)))
steps=[("1. Shell geometry fit to Migkas\n(A/A_all, D=58.2, R=14.3 Gpc)",G-1),
       ("2. + SR Doppler + central mass\n(β₁=0.587, q=0.230, μ=0.093)",R),
       ("3. Solve modes m=1–5\n(m=3–5 hold 92% of residual)",strip(R,5)),
       ("4. Solve modes m=1–6",strip(R,6)),
       ("5. Solve modes m=1–12",strip(R,12))]
r0=rms(G-1)
plt.rcParams.update({"font.size":15})
fig=plt.figure(figsize=(14,24)); gs=fig.add_gridspec(6,1,height_ratios=[1,1,1,1,1,1.3],hspace=0.95)
fig.suptitle("CMB flattening: model steps and error left",fontsize=24,y=0.93)
cols=["k","C3","C1","C2","C0"]
for i,(t,x) in enumerate(steps):
    ax=fig.add_subplot(gs[i]); r=rms(x); sc=10**np.floor(np.log10(np.abs(x).max()))
    ax.plot(ELL,x/sc,color=cols[i],lw=2.4)
    ax.axhline(0,color="gray",lw=.7); ax.set_xlim(0,360); ax.grid(alpha=.3)
    ax.set_ylabel(rf"$\Delta T/T\ \times10^{{{-int(np.log10(sc))}}}$")
    ax.set_title(t.replace("\n","  ")+(f"\nRMS {r:.2e}   |   starting point" if i==0 else f"\nRMS {r:.2e}   |   removed {100*(1-r/r0):.3f}% by amplitude, {100*(1-(r/r0)**2):.7f}% by variance"),fontsize=15)
    if i==4: ax.set_xlabel(r"Galactic longitude $\ell$ (deg)")
ax=fig.add_subplot(gs[5]); rs=[rms(x) for _,x in steps]
lab=["shell\nfit","+ SR +\nmass","m=1–5","m=1–6","m=1–12"]
ax.bar(range(5),rs,color=cols); ax.set_yscale("log"); ax.set_xticks(range(5)); ax.set_xticklabels(lab)
ax.axhline(1e-5,ls="--",color="purple",lw=2,label=r"CMB level $10^{-5}$")
for i,r in enumerate(rs): ax.text(i,r*1.4,f"{r:.1e}\n{r/1e-5:,.0f}× CMB" if r>1e-5 else f"{r:.1e}\n{r/1e-5:.2f}× CMB",ha="center",fontsize=13)
ax.set_ylim(1e-7,1); ax.set_ylabel("RMS anisotropy"); ax.grid(alpha=.3,axis="y")
ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),frameon=False); ax.set_title("Progress (log scale)",pad=40)
fig.savefig("cmb-progress.png",dpi=90,bbox_inches="tight")
print([f"{r:.3e}" for r in rs])

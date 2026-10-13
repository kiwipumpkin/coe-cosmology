import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
r=[5.751e-02,2.373e-03,6.113e-04,2.168e-04,1.833e-06]; cmb=1e-5
lab=["Shell fit\n(start)","+ SR +\nmass","m=1–5\nremoved","m=1–6\nremoved","m=1–12\nremoved"]
x=[v/cmb for v in r]; plt.rcParams.update({"font.size":15})
fig,ax=plt.subplots(2,1,figsize=(13,12),gridspec_kw=dict(hspace=1.25,height_ratios=[3,1]))
a=ax[0]; a.bar(range(5),x,color=["k","C3","C1","C2","C0"]); a.set_yscale("log"); a.set_ylim(0.05,3e4)
a.axhline(1,ls="--",color="purple",lw=2,label="CMB level (1×)")
for i,v in enumerate(x): a.text(i,v*1.3,f"{v:,.0f}×" if v>=1 else f"{v:.2f}×",ha="center",fontsize=15)
for i in range(1,5):
    a.text(i-0.5,1.2e4,f"÷{x[i-1]/x[i]:.1f}",ha="center",fontsize=15,color="dimgray",fontweight="bold")
a.set_xticks(range(5)); a.set_xticklabels(lab); a.set_ylabel("RMS ΔT/T in CMB units (log)")
a.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),frameon=False); a.set_title("Anisotropy left after each step",pad=40)
b=ax[1]; tot=x[0]-x[4]; parts=[x[i-1]-x[i] for i in range(1,5)]
names=["SR + mass","m=1–5","m=6","m=7–12"]; cols=["C3","C1","C2","C0"]; left=0
for p,n,cc in zip(parts,names,cols):
    b.barh(0,p,left=left,color=cc,label=f"{n}: {p:,.0f}× removed ({100*p/x[0]:.3f}% of start)"); left+=p
b.barh(0,x[4],left=left,color="purple",label=f"left: {x[4]:.2f}× ({100*x[4]/x[0]:.5f}%)")
b.set_xlim(0,x[0]*1.02); b.set_yticks([]); b.set_xlabel("share of the starting 5,751× (linear)")
b.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=1,frameon=False,fontsize=14); b.set_title("Where the starting anisotropy went (linear share)",pad=175)
fig.savefig("cmb-waterfall.png",dpi=95,bbox_inches="tight")

"""Brief (Oct 11): comparison charts for the mass-rising-outward shells, in the style of shells_remake/shells_compare_density.png
and structure_compare.png. Reads saved results only: /workspace/clump_scan/shells_remake/shells_massout.json."""
import json,numpy as np,matplotlib,matplotlib.ticker
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":12})
M=json.load(open("/workspace/clump_scan/shells_remake/shells_massout.json"))
N=12; R0=3.0; W0=1.0
def edges(q): w=W0*q**np.arange(N); return R0+np.r_[0,np.cumsum(w)]
cols={"constant":"#555555","widening":"#1747a6","narrowing":"#b03050"}; Rs={"constant":15.00,"widening":24.00,"narrowing":12.37}
fig,ax=plt.subplots(figsize=(10,6.5))
for k in ["constant","widening","narrowing"]:
    W=M["waves"][k]; e=edges(W["q"]); xe=e/e[-1]; rho=np.array(W["rho"]); r0=np.array(W["rho_equal_mass"])
    ax.step(np.r_[xe[:-1],xe[-1]],np.r_[rho,rho[-1]],where="post",color=cols[k],lw=2.4,label=f"{k}, mass rising outward (R = {Rs[k]:.2f})")
    ax.step(np.r_[xe[:-1],xe[-1]],np.r_[r0,r0[-1]],where="post",color=cols[k],lw=1.2,ls=":",label=f"{k}, equal mass")
ax.set_yscale("log"); ax.set_xlim(0,1); ax.set_xlabel("normalised radius  r / R"); ax.set_ylabel("shell density ÷ innermost equal-mass shell"); ax.grid(alpha=0.3,which="both")
ax.set_title("Mass rising outward (3×): density profile for the three spacings",fontsize=14,fontweight="bold")
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.14),ncol=2,frameon=False,fontsize=10.5)
fig.savefig("figures/new/shells_massout_compare_density.png",dpi=130,bbox_inches="tight"); plt.close(fig)
fig,ax=plt.subplots(2,1,figsize=(10,9),sharex=True); fig.subplots_adjust(hspace=0.12)
for k in ["constant","widening","narrowing"]:
    S=M["structure"][k]; e=edges(M["waves"][k]["q"]); xe=e/e[-1]
    s=[r["sigma_measured"] for r in S["rows"]]; sp=[r["spacing_measured"] or np.nan for r in S["rows"]]
    ax[0].step(np.r_[xe[:-1],xe[-1]],np.r_[s,s[-1]],where="post",color=cols[k],lw=2.4,label=f"{k} (R = {Rs[k]:.2f})")
    ax[1].step(np.r_[xe[:-1],xe[-1]],np.r_[sp,sp[-1]],where="post",color=cols[k],lw=2.4)
ax[0].set_ylabel("clumpiness rms(ρ/ρ_shell − 1)"); ax[0].set_ylim(0,1.6); ax[0].grid(alpha=0.3)
ax[1].set_ylabel("typical clump spacing (radius units)"); ax[1].set_yscale("log"); ax[1].set_yticks([0.3,0.5,1,2]); ax[1].set_yticklabels(["0.3","0.5","1","2"]); ax[1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[1].set_ylim(0.3,2.3); ax[1].grid(alpha=0.3,which="both"); ax[1].set_xlim(0,1); ax[1].set_xlabel("normalised radius  r / R")
ax[0].set_title("Mass rising outward (3×): clumpiness and clump spacing vs radius",fontsize=14,fontweight="bold")
fig.legend(*ax[0].get_legend_handles_labels(),loc="lower center",ncol=3,frameon=False,fontsize=11,bbox_to_anchor=(0.5,0.0)); fig.subplots_adjust(bottom=0.13)
fig.savefig("figures/new/structure_massout_compare.png",dpi=130,bbox_inches="tight"); plt.close(fig)
print({k:[round(r["sigma_measured"],2) for r in M["structure"][k]["rows"]][-1] for k in cols})

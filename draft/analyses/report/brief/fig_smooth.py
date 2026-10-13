"""Ring-plane ('side view', Chris's convention) panel of the 1x smooth field at 2 Gpc smoothing. Reads saved results only."""
import os,sys,json,numpy as np,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
OUT=os.path.abspath("figures"); SF="/workspace/clump_scan/inside_visible/smooth_field"
sys.path.insert(0,SF); sys.path.insert(0,os.path.dirname(SF)); os.chdir(SF); sys.argv=[""]
from sf_core import Grid
from inside_core import C_G
plt.rcParams.update({"font.size":8,"axes.titlesize":8.5,"savefig.dpi":300,"pdf.fonttype":42})
T=json.load(open("tradeoff.json"))
g=Grid(0.3); d=np.load("delta_H0.3_s2.npy"); F=np.full(g.mask.shape,np.nan); F[g.mask]=d
i0=np.argmin(np.abs(g.ax)); ext=[g.ax[0],g.ax[-1],g.ax[0],g.ax[-1]]; vm=0.05
fig,ax=plt.subplots(figsize=(3.4,3.0))
im=ax.imshow(F[:,:,i0].T,origin="lower",extent=ext,cmap="RdBu_r",vmin=-vm,vmax=vm)
ax.add_patch(Circle((0,0),14.3,fill=False,ls="--",color="k",lw=1)); ax.plot(0,0,"k+",ms=8,mew=1.5)
u=C_G/np.linalg.norm(C_G); ax.annotate("",xy=(u[0]*15.3,u[1]*15.3),xytext=(u[0]*12.0,u[1]*12.0),arrowprops=dict(arrowstyle="->",color="0.2",lw=1))
ax.text(u[0]*16.3,u[1]*16.3-2.3,"to CoE",fontsize=7,color="0.2",ha="center")
ax.set_xlim(-19.5,15.5); ax.set_ylim(-15.5,15.5); ax.set_aspect("equal"); ax.set_xlabel("x (Gpc)"); ax.set_ylabel("y (Gpc)")
ax.set_title("side view: slice through the ring plane\nrms $\\delta$ = %.3f, peak $|\\delta|$ = %.3f"%(T["H0.3_s2"]["at1x"]["rms"],T["H0.3_s2"]["at1x"]["peak"]))
cb=fig.colorbar(im,ax=ax,shrink=0.9,pad=0.03); cb.set_label("$\\delta$ (blue under-, red over-dense)",fontsize=7); cb.ax.tick_params(labelsize=6.5)
fig.savefig(f"{OUT}/fig8a_smooth_field_side.pdf",bbox_inches="tight"); fig.savefig(f"{OUT}/fig8a_smooth_field_side.png",bbox_inches="tight",dpi=200)
print("ok")

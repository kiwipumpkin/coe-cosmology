import json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from lcdm import sigma_G, sigma_TH
plt.rcParams.update({"font.size":14,"axes.titlesize":15,"legend.fontsize":12})
T=json.load(open("tradeoff.json"))
COL={0:"#555555",1:"#c0392b",2:"#8e44ad",4:"#2471a3"}
LAB={0:"no extra smoothing (0.3 Gpc grid)",1:"smoothed 1 Gpc",2:"smoothed 2 Gpc",4:"smoothed 4 Gpc"}
fig,ax=plt.subplots(figsize=(10.5,7.5))
for s in [0,1,2,4]:
    c=T[f"H0.3_s{s}"]["curve"]; x=np.array([q["rms"] for q in c]); y=np.array([q["left"] for q in c]); o=np.argsort(x)
    ax.plot(x[o],y[o],"-",color=COL[s],lw=2.6,label=LAB[s])
    a=T[f"H0.3_s{s}"]["at1x"]; ax.plot(a["rms"],1.0,"o",color=COL[s],ms=9,zorder=5)
    R=0.3 if s==0 else s; lo,hi=sigma_G(R),sigma_TH(R)
    ax.axvspan(lo,hi,color=COL[s],alpha=0.18)
ax.axhline(1.0,color="green",lw=2); ax.text(1.3e-5,1.12,"CMB flat target 1×",color="green",fontsize=12)
ax.axhline(237.3,color="0.5",ls="--",lw=1); ax.text(1.3e-5,262,"start 237×",color="0.4",fontsize=11)
ax.text(1.3e-5,0.012,"shaded strips = what ΛCDM expects (rms δ today)\nat the same smoothing scale, same colours",fontsize=11,color="0.3")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-5,0.3); ax.set_ylim(0.008,500)
ax.set_xlabel("rms density contrast δ inside the sphere"); ax.set_ylabel("leftover (× 10⁻⁵ floor)")
ax.set_title("Smooth density field inside our visible sphere: how strong must it be?")
ax.legend(loc="lower left",bbox_to_anchor=(0.0,0.12),frameon=False); fig.tight_layout(); fig.savefig("tradeoff_leftover_vs_rms_delta.png",dpi=110); plt.close(fig)

# ---- map of the 1x field, smoothing 2 Gpc
import sys; sys.argv=[""]
from sf_core import Grid
from inside_core import C_G
g=Grid(0.3); s=2; d=np.load(f"delta_H0.3_s{s}.npy"); F=np.full(g.mask.shape,np.nan); F[g.mask]=d
i0=np.argmin(np.abs(g.ax)); ext=[g.ax[0],g.ax[-1],g.ax[0],g.ax[-1]]; vm=0.05
fig,axs=plt.subplots(1,2,figsize=(15.5,7.4))
for ax,(img,xl,yl,tt) in zip(axs,[(F[:,:,i0].T,"x (Gpc)","y (Gpc)","top view: slice through the ring plane"),
                                  (F[:,i0,:].T,"x (Gpc)","z (Gpc)","side view: vertical slice along ℓ = 0°–180° (CoE is at 195°)")]):
    im=ax.imshow(img,origin="lower",extent=ext,cmap="RdBu_r",vmin=-vm,vmax=vm)
    ax.add_patch(Circle((0,0),14.3,fill=False,ls="--",color="k",lw=1.5)); ax.plot(0,0,"k+",ms=16,mew=2.5)
    ax.set_xlim(-17,15.5); ax.set_ylim(-15.5,15.5); ax.set_aspect("equal"); ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(tt,fontsize=13)
u=C_G/np.linalg.norm(C_G); axs[0].annotate("",xy=(u[0]*15.3,u[1]*15.3),xytext=(u[0]*12.0,u[1]*12.0),arrowprops=dict(arrowstyle="->",color="0.2",lw=1.5))
axs[0].text(u[0]*15.0,u[1]*15.0-1.3,"to CoE",fontsize=11,color="0.2",ha="center")
cb=fig.colorbar(im,ax=axs,shrink=0.85,pad=0.02); cb.set_label("δ   (blue = under-dense, red = over-dense)")
fig.suptitle("Weakest smooth field (2 Gpc smoothing) that reaches 1×: rms δ = %.3f, peak |δ| = %.3f"%(T["H0.3_s2"]["at1x"]["rms"],T["H0.3_s2"]["at1x"]["peak"]),fontsize=14)
fig.savefig("map_minimal_field_2Gpc.png",dpi=110,bbox_inches="tight"); plt.close(fig)

# ---- sky direction: delta vs longitude just inside the ring, in the plane
P=g.P; r=np.linalg.norm(P,axis=1); lon=np.degrees(np.arctan2(P[:,1],P[:,0]))%360; sel=(r>11)&(np.abs(P[:,2])<1.5)
b=np.arange(0,361,5); cen=b[:-1]+2.5
fig,ax=plt.subplots(figsize=(11,5.8))
prof=[d[sel&(lon>=a)&(lon<a+5)].mean() for a in b[:-1]]
ax.plot(cen,np.array(prof),"-",color="#8e44ad",lw=2.8)
ax.axhline(0,color="0.5",lw=1)
for L in (75,162,228,315): ax.axvline(L,color="#c0392b",ls=":",lw=2)
ax.text(75,0.047,"old four-mass\nlongitudes (outside)",color="#c0392b",fontsize=11,ha="center",va="top")
ax.set_xlim(0,360); ax.set_ylim(-0.05,0.05); ax.set_xticks(range(0,361,45))
ax.set_xlabel("sky longitude ℓ (deg)"); ax.set_ylabel("δ (mean, 11–13.8 Gpc, near ring plane)")
ax.set_title("Where the 1× field is dense or empty, just inside the CMB ring (2 Gpc smoothing)")
fig.tight_layout(); fig.savefig("sky_delta_vs_longitude.png",dpi=110); plt.close(fig)

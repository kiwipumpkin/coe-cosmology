import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
from matplotlib.colors import ListedColormap
import sys; sys.argv=[sys.argv[0]]
from cb_core import Q_SG, LAM, AXIS
RV=14.3; cm=ListedColormap(["#3b6fd8","#d62728"])
def chess(x,y):
    P=np.c_[x.ravel(),y.ravel(),np.zeros(x.size)]@Q_SG; c=LAM/2
    return ((np.floor(P[:,0]/c)+np.floor(P[:,1]/c)+np.floor(P[:,2]/c))%2==0).reshape(x.shape).astype(float)
def img(ax,cx,cy,half,n):
    g=np.linspace(-half,half,n); X,Y=np.meshgrid(cx+g,cy+g); S=chess(X,Y)
    ax.imshow(S,origin="lower",extent=(cx-half,cx+half,cy-half,cy+half),cmap=cm,interpolation="nearest",vmin=0,vmax=1)
a=np.radians(195.16); pc=np.array([np.cos(a),np.sin(a)])*RV
fig=plt.figure(figsize=(11,19)); gs=fig.add_gridspec(3,3,height_ratios=[2.6,1,0.25],hspace=0.12,wspace=0.12)
ax=fig.add_subplot(gs[0,:]); img(ax,0,0,15.8,3800)
mask=np.hypot(*np.meshgrid(np.linspace(-15.8,15.8,600),np.linspace(-15.8,15.8,600)))>RV
ax.imshow(np.where(mask,1.0,np.nan),origin="lower",extent=(-15.8,15.8,-15.8,15.8),cmap=ListedColormap(["white"]),alpha=0.55)
ax.add_patch(Circle((0,0),RV,fill=False,color="#1a9e3f",lw=2.5))
ax.plot(0,0,"k*",ms=16,mec="w"); ax.text(0.3,0.4,"us",fontsize=12,color="k",weight="bold",bbox=dict(fc="w",alpha=0.8,lw=0))
ax.annotate("",xy=(pc*1.03),xytext=pc*0.55,arrowprops=dict(arrowstyle="-|>",color="k",lw=2.5))
ax.text(*(pc*0.66+np.array([0,0.7])),"to CoE (l = 195°, 58.2 Gpc)",fontsize=11,ha="center",bbox=dict(fc="w",alpha=0.85,lw=0))
ax.annotate("",xy=(RV,-0.6),xytext=(0,-0.6),arrowprops=dict(arrowstyle="<->",color="k",lw=1.5))
ax.text(RV/2,-1.3,"14.3 Gpc = 167 chess cells (83 full periods) to the ring",ha="center",fontsize=11,bbox=dict(fc="w",alpha=0.85,lw=0))
ax.text(0,RV+0.35,"334 cells across our visible sphere",ha="center",fontsize=11,bbox=dict(fc="w",alpha=0.85,lw=0))
sizes=[3.0,1.0,0.3]; ax.add_patch(Rectangle(pc-sizes[0]/2,sizes[0],sizes[0],fill=False,color="k",lw=2))
ax.set_xlim(-15.8,15.8); ax.set_ylim(-15.8,15.8); ax.set_aspect("equal"); ax.set_xlabel("Gpc (ring plane = galactic plane, x toward l = 0°)"); ax.set_ylabel("Gpc")
ax.set_title("Side view: Einasto-style 3D chessboard through our visible universe (slice through the ring plane)\ncells 85.7 Mpc = ½ of the 171 Mpc (120 h⁻¹ Mpc) period; grid axes = supergalactic axes as in Einasto Fig 1.\nAt this scale the 334 cells blur together; zooms below are centred on the ring point toward the CoE",fontsize=11)
for i,s in enumerate(sizes):
    ax2=fig.add_subplot(gs[1,i]); img(ax2,pc[0],pc[1],s/2,900)
    t=np.linspace(0,2*np.pi,20000); ax2.plot(RV*np.cos(t),RV*np.sin(t),color="#1a9e3f",lw=3)
    if i<2: ax2.add_patch(Rectangle(pc-sizes[i+1]/2,sizes[i+1],sizes[i+1],fill=False,color="k",lw=2))
    ax2.set_xlim(pc[0]-s/2,pc[0]+s/2); ax2.set_ylim(pc[1]-s/2,pc[1]+s/2); ax2.set_aspect("equal")
    lab={3.0:"3 Gpc box (~35 cells across)",1.0:"1 Gpc box (~12 cells)",0.3:"300 Mpc box (~3.5 cells)"}[s]
    ax2.set_title(f"zoom {i+1}: {lab}",fontsize=10); ax2.tick_params(labelsize=8)
    ax2.annotate("",xy=pc+(pc/RV)*s*0.45,xytext=pc+(pc/RV)*s*0.15,arrowprops=dict(arrowstyle="-|>",color="k",lw=2))
    ax2.text(*(pc+(pc/RV)*s*0.3+np.array([0,s*0.06])),"CoE",fontsize=9,ha="center",bbox=dict(fc="w",alpha=0.8,lw=0))
    sb={3.0:1.0,1.0:0.2,0.3:0.0857}[s]; x0=pc[0]-s*0.45; y0=pc[1]-s*0.43
    ax2.plot([x0,x0+sb],[y0,y0],"k-",lw=4); ax2.text(x0+sb/2,y0+s*0.03,f"{sb*1000:.0f} Mpc" if sb<1 else "1 Gpc",ha="center",fontsize=9,bbox=dict(fc="w",alpha=0.85,lw=0))
axl=fig.add_subplot(gs[2,:]); axl.axis("off")
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
axl.legend([Patch(color="#d62728"),Patch(color="#3b6fd8"),Line2D([],[],color="#1a9e3f",lw=3),Line2D([],[],marker="*",color="k",lw=0,ms=12),Patch(fc="w",ec="k")],
  ["over-dense cell (+δ)","under-dense cell (−δ)","CMB ring = edge of our 14.3 Gpc sphere","us","next zoom box"],loc="center",ncol=3,fontsize=10,frameon=False)
fig.text(0.5,0.105,"The ring plane is tilted 6.3° from a face of the supergalactic grid, so the slice cuts through neighbouring layers: cells appear as squares turned ~47° (diamonds, because the supergalactic X axis\npoints to l = 137°) whose layout slowly shifts by one cell every ≈0.8 Gpc (the faint bands in zoom 1). Outside our sphere (faded) the same pattern continues.",ha="center",fontsize=9)
fig.savefig("fig3_side_view_checkerboard.png",dpi=150,bbox_inches="tight")

import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Patch
from matplotlib.lines import Line2D
import sys; sys.argv=[sys.argv[0]]
from ri_core import LAM, C_G
RV=14.3; cmap=plt.get_cmap("RdBu_r"); C2=C_G[:2]
def img(ax,cx,cy,half,n,lam=LAM):
    g=np.linspace(-half,half,n); X,Y=np.meshgrid(cx+g,cy+g)
    r=np.hypot(X-C2[0],Y-C2[1]); S=np.cos(2*np.pi*r/lam); S[np.hypot(X,Y)>RV]=np.nan
    cm=cmap.copy(); cm.set_bad("white")
    ax.imshow(S,origin="lower",extent=(cx-half,cx+half,cy-half,cy+half),cmap=cm,vmin=-1,vmax=1,interpolation="nearest")
a=np.radians(195.16); pc=np.array([np.cos(a),np.sin(a)])*RV
fig=plt.figure(figsize=(11,19)); gs=fig.add_gridspec(3,3,height_ratios=[2.6,1,0.25],hspace=0.12,wspace=0.12)
ax=fig.add_subplot(gs[0,:]); img(ax,0,0,15.8,3800)
ax.add_patch(Circle((0,0),RV,fill=False,color="#1a9e3f",lw=2.5))
ax.plot(0,0,"k*",ms=16,mec="w"); ax.text(0.3,0.4,"us",fontsize=12,weight="bold",bbox=dict(fc="w",alpha=0.8,lw=0))
ax.annotate("",xy=pc*1.03,xytext=pc*0.55,arrowprops=dict(arrowstyle="-|>",color="k",lw=2.5))
ax.text(*(pc*0.66+np.array([0,0.7])),"to CoE (l = 195°, 58.2 Gpc)",fontsize=11,ha="center",bbox=dict(fc="w",alpha=0.85,lw=0))
u=-pc/RV; ax.annotate("",xy=-pc,xytext=pc,arrowprops=dict(arrowstyle="<->",color="k",lw=1.2))
ax.text(4,-2.2,"along the CoE axis: 43.9 → 72.5 Gpc from the CoE\n= 28.6 Gpc / 171 Mpc ≈ 167 crests across our sphere",fontsize=10.5,ha="center",bbox=dict(fc="w",alpha=0.85,lw=0))
ax.text(0,RV+0.4,"perpendicular to the axis the shells curve away: only ≈10 crests from centre to edge (1.7 Gpc of CoE distance)",fontsize=9.5,ha="center",bbox=dict(fc="w",alpha=0.85,lw=0))
sizes=[3.0,1.0,0.3]; ctr=lambda s_: pc*(1-0.3*s_/RV); ax.add_patch(Rectangle(ctr(sizes[0])-sizes[0]/2,sizes[0],sizes[0],fill=False,color="k",lw=2))
ax.set_xlim(-15.8,15.8); ax.set_ylim(-15.8,15.8); ax.set_aspect("equal"); ax.set_xlabel("Gpc (ring plane = galactic plane, x toward l = 0°)"); ax.set_ylabel("Gpc")
ax.set_title("Side view: radial blast-wave shells around the CoE, drawn only inside our visible sphere\nfull wavelength 171 Mpc (Einasto 120 h⁻¹ Mpc, h = 0.7); red over-dense, blue under-dense.\nAt this scale the 167 crests blur together; zooms below are centred on the ring point toward the CoE",fontsize=11)
for i,s in enumerate(sizes):
    ax2=fig.add_subplot(gs[1,i]); q=ctr(s); img(ax2,q[0],q[1],s/2,900)
    t=np.linspace(0,2*np.pi,20000); ax2.plot(RV*np.cos(t),RV*np.sin(t),color="#1a9e3f",lw=3)
    if i<2: ax2.add_patch(Rectangle(ctr(sizes[i+1])-sizes[i+1]/2,sizes[i+1],sizes[i+1],fill=False,color="k",lw=2))
    ax2.set_xlim(q[0]-s/2,q[0]+s/2); ax2.set_ylim(q[1]-s/2,q[1]+s/2); ax2.set_aspect("equal")
    lab={3.0:"3 Gpc box (~17 crests)",1.0:"1 Gpc box (~6 crests)",0.3:"300 Mpc box (~1.75 waves)"}[s]
    ax2.set_title(f"zoom {i+1}: {lab}",fontsize=10); ax2.tick_params(labelsize=8)
    ax2.annotate("",xy=q+(pc/RV)*s*0.25,xytext=q-(pc/RV)*s*0.1,arrowprops=dict(arrowstyle="-|>",color="k",lw=2))
    ax2.text(*(q+(pc/RV)*s*0.1+np.array([0,s*0.08])),"to CoE",fontsize=9,ha="center",bbox=dict(fc="w",alpha=0.8,lw=0))
    sb={3.0:1.0,1.0:0.2,0.3:0.1714}[s]; x0=q[0]+s*0.05; y0=q[1]-s*0.43
    ax2.plot([x0,x0+sb],[y0,y0],"k-",lw=4); ax2.text(x0+sb/2,y0+s*0.03,{1.0:"1 Gpc",0.2:"200 Mpc",0.1714:"171 Mpc = 1 wave"}[sb],ha="center",fontsize=9,bbox=dict(fc="w",alpha=0.85,lw=0))
axl=fig.add_subplot(gs[2,:]); axl.axis("off")
axl.legend([Patch(color=cmap(0.95)),Patch(color=cmap(0.05)),Line2D([],[],color="#1a9e3f",lw=3),Line2D([],[],marker="*",color="k",lw=0,ms=12),Patch(fc="w",ec="k"),Patch(fc="w",ec="0.7")],
  ["over-dense shell (+δ)","under-dense shell (−δ)","CMB ring = edge of our 14.3 Gpc sphere","us","next zoom box","outside our sphere: no wave in this model"],loc="center",ncol=3,fontsize=10,frameon=False)
fig.text(0.5,0.105,"Shell spacing is the same everywhere (171 Mpc in CoE distance); the shells are spheres 44–73 Gpc in radius, so inside our sphere they look like nearly parallel,\ngently curved bands crossing the sightline to the CoE. Near l = 195° and 15° the ring runs along the shells (it stays in one shell for a while); near l = 105° and 285° it cuts straight across them.",ha="center",fontsize=9)
fig.savefig("fig1_side_view_radial_shells.png",dpi=150,bbox_inches="tight")

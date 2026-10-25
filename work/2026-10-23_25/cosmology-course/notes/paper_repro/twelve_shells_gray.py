import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, Rectangle
c=1.0; Rrec=3.0; D=lambda t: 3.0+1.0*t     # D1=4,D2=5,D3=6
import matplotlib.colors as mcolors
import sys; K=int(sys.argv[1]); P=np.load('p_r2.npy'); p,rho0=P[2*(K-1)],P[2*(K-1)+1]; M=lambda r: r**p/rho0
allrho=[]
for tt in (1,2,3):
    L=np.linspace(D(tt)+c*tt,D(tt)-c*tt,2*K*tt+1)
    allrho+= [M((L[i-1]+L[i])/2)/((4/3)*np.pi*(L[i-1]**3-L[i]**3)) for i in range(1,2*K*tt+1)]
rmax=1.0; norm=mcolors.Normalize(0.95,1.0); cmap=plt.cm.Greys
fig,axs=plt.subplots(1,3,figsize=(24,10))
for ax,t in zip(axs,[1,2,3]):
    d=D(t); R=c*t; n=2*K*t
    lines=np.linspace(d+R,d-R,n+1)          # n+1 boundary lines -> n bands
    assert abs(lines[-1]-(d-R))<1e-9 and abs(lines[0]-(d+R))<1e-9
    assert abs((lines[0]-lines[-1])-2*R)<1e-9 and abs(lines[-1]-Rrec)<1e-12 and abs(Rrec-(d-R))<1e-12 and len(lines)==n+1
    if t>1: assert D(t-1)<d
    for i in range(1,n+1):
        ro,ri=lines[i-1],lines[i]; m=(ro+ri)/2
        ax.add_patch(Wedge((0,0),ro,-90,90,width=ro-ri,fc=cmap(norm(M(m)/((4/3)*np.pi*(ro**3-ri**3))/rmax)),ec="none",alpha=0.55,zorder=0))
        pass #,ri,ro,f"{1/((4/3)*np.pi*(ro**3-ri**3)):.5f}",f"{1/((4/3)*np.pi*(ro**3-ri**3))/rmax:.3f}")
        ax.text(-0.25,-m,str(i),ha="right",va="center",fontsize=13 if K==1 else 8,weight="bold")
        for ang in np.deg2rad([-60,-30,0,30,60]):
            u=np.array([np.cos(ang),np.sin(ang)]); ln=0.08*m
            ax.annotate('',tuple((m+ln/2)*u),tuple((m-ln/2)*u),arrowprops=dict(arrowstyle='->',color='navy',lw=1.7),zorder=2)

    for j,s_ in enumerate(lines):
        ax.add_patch(Wedge((0,0),s_,-90,90,width=1e-4,ec="k",lw=2 if j in(0,n) else 1.0))
    ax.add_patch(Wedge((0,0),Rrec,-90,90,fc="white",ec="k",hatch="////",lw=1.5))
    ax.text(1.2,-1.6,"R_rec = D−R",va="center",ha="center",bbox=dict(fc='white',ec='none'),fontsize=12,weight="bold"); ax.text(-0.25,-0.35,"CoE",ha="right",fontsize=12)
    centres=[]
    for k in range(t-1,0,-1):                       # concentric on current OBS
        pass
    pass
    if t>1: ax.text(d,-(t-1)+0.3,"already seen",ha="center",fontsize=9,weight="bold",bbox=dict(fc="white",alpha=0.8,ec="none"))
    ROY=["red","orange","gold","green","blue","indigo","violet"]
    for k in range(t,0,-1):
        ax.add_patch(Circle((d,0),c*k,fc=ROY[k-1],alpha=0.3,ec="none",zorder=3)); centres.append(d)
        ax.add_patch(Circle((d,0),c*k,fc="none",ec=ROY[k-1],lw=2.5,zorder=3))
    assert len(centres)==t and all(abs(x-d)<1e-12 for x in centres)
    ax.text(d,R+0.15,"seeing now (radius R = c·t)",ha="center",color="k",fontsize=12,bbox=dict(fc="white",ec="none",alpha=0.85))
    for x in (d-R,d+R): ax.plot(x,0,"o",color="tab:blue",ms=6)
    ax.annotate("",(d+R,-R-0.35),(d-R,-R-0.35),arrowprops=dict(arrowstyle="<->",color="tab:blue"))
    ax.text(d,-R-0.85,"2R",ha="center",color="tab:blue",fontsize=13,weight="bold")
    ax.plot(d,0,"ko",ms=9); ax.text(d+0.1,0.25,"OBS",fontsize=13,weight="bold")
    ax.annotate("",(d+0.8,0),(d,0),arrowprops=dict(arrowstyle="->",lw=2.5))
    ax.plot([0,11],[0,0],color="gray",lw=0.8); ax.axvline(0,color="k",lw=2)
    y=-(d+R)-0.8; ax.annotate("",(d,y),(0,y),arrowprops=dict(arrowstyle="<->")); ax.plot([d,d],[y,0],"k:",lw=1)
    ax.text(d/2,y-0.5,f"D_{t} = {d:.0f}",ha="center",fontsize=13)
    ax.set_xlim(-2,11); ax.set_ylim(-11,10); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(f"t = {t}   ({n} shells)",fontsize=17)
fig.suptitle("Round shells centred on the CoE. The circle we see around OBS touches the youngest shell (D−R) on the CoE side and shell 1 (D+R) on the far side,\n"
 "so the shells span 2R across our circle. Each layer emitted once at 3000 K crossing R_rec. Shells are the bands between lines (n+1 lines, n bands). Navy arrows: shell speed (outer faster). Same scale, units of c·Δt.",fontsize=13)
fig.text(0.5,0.02,f"shell 1: oldest, coolest, fastest     shell {6*K}: youngest, warmest",ha="center",fontsize=16,weight="bold")
plt.subplots_adjust(left=0.02,right=0.93,top=0.87,bottom=0.06,wspace=0.05)
cax=fig.add_axes([0.945,0.15,0.012,0.65]); cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_label(f"density ÷ innermost band (content ∝ r^{p:.3f} ÷ shell volume); shared scale",fontsize=12)
plt.savefig("twelve_shells_gray.png",dpi=200)
